module CellularAutomata
using Distributions
using Random
using StatsBase

export colonise,coloniseSuitWeight,coloniseNeighbourWeight,extinction,extinctionNeighbourWeight,OccurenceCellularAutomata,OccurrenceCellularAutomataNB

# struct DispersalParams
#     meanDispersal::Float64
#     dispersalProbability::Float64
#     maxNumberDispersers::Int64
# end
# struct NeighbourhoodParams
#     neighbourSurvivalWeight::Float64
#     neighbourDispersalWeight::Float64
#     neighbourWeightMatrix::Matrix{Float64}
#     r::Int64
# end
mutable struct SpeciesEnvironment
    suitabilityInitial::Matrix{Float64}
    suitabilityActive::Matrix{Float64}
end
mutable struct SpeciesCellularAutomataNeighbourEffect
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    # Dispersal parameters
    meanDispersal::Float64
    dispersalProbability::Float64
    maxNumberDispersers::Int64
    # Weight parameters
    neighbourSurvivalWeight::Float64
    neighbourDispersalWeight::Float64
    neighbourWeightMatrix::Matrix{Float64}
    r::Int64
end
mutable struct SpeciesCellularAutomata
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    # Dispersal parameters
    meanDispersal::Float64
    dispersalProbability::Float64
    maxNumberDispersers::Int64
end

"""
    newPos(meanDistance,cartIdx)

Calculates a new position that does not include a starting position (cartIdx).
This is achived stochastically by sampling distance from an exponential
distribution and direction of movement from a uniform distribution between 0-360.

Note: The selected distance as a constant added to it in order to ensure a
      distance less than 1 does not select a position the same as the starting one.
      This value is +0.7071 because this is the ratio between the distance to
      a direct:diagonal neighbour.
      This ensures that rounding the X and Y values the direct neighbours are
      selected more frequently than the diagonals. Higher values will bias to
      diagonals and lower values will result in selection of the starting position.

"""
function newPos(meanDistance::Float64,startPos::CartesianIndex,rng)
    distance = rand(rng,Exponential(meanDistance),1)
    # Add constant to ensure selection of pos outside current cell.
    distance = distance[1]+0.7071
    angle = 360.0*rand()
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + startPos[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + startPos[1]
    coord = (y,x)
    return(coord)
end

"""
    selectProportion(pa,caIndex,dispersalProbability)

Given a 2D array of occurence, randomly select cells that will disperse at
the next time step. The number of dispersers is determined stochastically
by sampling from a binomial distribution.

"""
function selectProportion(pa::Matrix{Float64},caIndex::Array{CartesianIndex{2},2},probability::Float64)
    nPresences = Int(sum(pa))
    total = Array{CartesianIndex}(undef, nPresences)
    counter = 1
    for idx in caIndex
        if pa[idx] === 1.0
            total[counter] = idx
            counter+=1
        end

    end
    dispersers = sample(total,rand(Binomial(nPresences,probability)))
end
"""
    neighbourHoodWeight(ca)

Select a single cell to colonise

"""
function neighbourHoodWeight(pa::Matrix{Float64},suit::Matrix{Float64},idx::CartesianIndex,shape::Tuple,weightMatrix::Matrix{Float64},r::Int64)
    # Assumes square neighbourhood
    ymin = idx[1]-r
    ymax = idx[1]+r
    xmin = idx[2]-r
    xmax = idx[2]+r
    idxYMin = 1
    idxYMax = 1+2*r
    idxXMin = 1
    idxXMax = 1+2*r
    # Check boundary conditions
    if ymin < 1
        ymin = 1
        idxYMin = r+(2-idx[1])
    end
    if ymax > shape[1]
        ymax = shape[1]
        idxYMax = (r+1)+shape[1]-idx[1]
    end
    if xmin < 1
         xmin = 1
         idxXMin = r+(2-idx[2])
      end
    if xmax > shape[2]
        xmax = shape[2]
        idxXMax = (r+1)+shape[2]-idx[2]
    end
    oc = vec(pa[ymin:ymax,xmin:xmax])
    su = vec(suit[ymin:ymax,xmin:xmax])
    weight = oc.*su.*vec(weightMatrix[idxYMin:idxYMax,idxXMin:idxXMax])
    #weight = oc.*vec(weightMatrix[idxYMin:idxYMax,idxXMin:idxXMax])
    return(sum(weight))
end
"""
    colonise(ca)

Select a single cell to colonise

"""
function colonise(ca::OccurrenceCellularAutomata)
    rng = MersenneTwister()
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        newXY = newPos(ca.meanDispersal,i,rng)
        if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
            ca.pa[newXY[1],newXY[2]] = 1.0
        end
    end
end

"""
    coloniseSuitWeight(ca)

Randomly selects a number cells that attempt to colonise nearby cells.
A maximum number of dispersers is specified and each cell selected to disperse
has


"""
function coloniseSuitWeight(ca::OccurrenceCellularAutomata,env::SpeciesEnvironment)
    rng = MersenneTwister()
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        # Determine maximum number of dispersers scaled by suitabiility
        meanNumDispersers = ca.maxNumberDispersers*env.suitabilityActive[i]
        if meanNumDispersers > 0.0
            numberDispersers = rand(rng,Poisson(meanNumDispersers))
            for j in 1:numberDispersers
                newXY = newPos(ca.meanDispersal,i,rng)
                if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
                    ca.pa[newXY[1],newXY[2]] = 1.0
                end
            end
        end
    end
end
"""
    coloniseNeighbourWeight(ca)

Randomly selects a number cells that attempt to colonise nearby cells.
Each cell selected to disperse has a random number of attempts at colonising
selected from a Poisson distribution. For each attempt, the new position is
selected randomly and converted to occupied (1.0).

Note: As this is neighbour weighted, the order of cells are randomised to prevent
      any bias due to the iterating over cells in the same order each time the
      function is called. This effectively makes the neighbourhood effects
      asynchronus across the landscape.
      An alternative would be to copy the occurence matrix and refer to this for
      neighbourhood weights. This would effectively make the neighbourhood effects
      synchronus (i.e. all processes happen at the same time, instantly).
"""
function coloniseNeighbourWeight(ca::OccurrenceCellularAutomataNB,env::SpeciesEnvironment)
    rng = MersenneTwister()
    shp = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    # Random shuffle to avoid grid index bias due to order of applying this function
    idxShuffle = sample(collect(1:1:length(dCells)),length(dCells),replace=false)
    for i in idxShuffle
        cellIndex = dCells[i]
        neighbourWeight = neighbourHoodWeight(ca.pa,env.suitabilityActive,cellIndex,shp,ca.neighbourWeightMatrix,ca.r)
        dispWeight = neighbourWeight * ca.neighbourDispersalWeight
        # Determine maximum number of dispersers scaled by suitabiility
        dispersalMultiplier = 1/(1+MathConstants.e^(-10*(env.suitabilityActive[cellIndex]-(1-dispWeight))))
        meanNumDispersers = (env.suitabilityActive[cellIndex] + (1-env.suitabilityActive[cellIndex])*dispersalMultiplier)*ca.maxNumberDispersers

        numberDispersers = rand(rng,Poisson(meanNumDispersers))
        for j in 1:numberDispersers
            newXY = newPos(ca.meanDispersal,cellIndex,rng)
            if newXY[2]>=1 && newXY[2] <= shp[2] && newXY[1] >=1 && newXY[1]<=shp[1]
                ca.pa[newXY[1],newXY[2]] = 1.0
            end
        end
    end
end


function extinction(ca::OccurrenceCellularAutomata)
    rng = MersenneTwister()
    for idx in ca.pa_cart_index
        if ca.pa[idx] === 1.0
            survived = rand(rng,Bernoulli(ca.suitability[idx]),1)
            if survived[1] == false
                ca.pa[idx] = 0
            end
        end
    end
end
"""
    extinctionNeighbourWeight(ca)

Each occupied cell is

Note: As this is neighbour weighted, the order of cells are randomised to prevent
      any bias due to the iterating over cells in the same order each time the
      function is called. This effectively makes the neighbourhood effects
      asynchronus across the landscape.
      An alternative would be to copy the occurence matrix and refer to this for
      neighbourhood weights. This would effectively make the neighbourhood effects
      synchronus (i.e. all processes happen at the same time, instantly).
"""
function extinctionNeighbourWeight(ca::OccurrenceCellularAutomataNB,env::SpeciesEnvironment)
    rng = MersenneTwister()
    # Randomise cell index to prevent bias from the neighbourhood weighting
    # Alternative could be to copy the the pa array.
    shp = size(ca.pa)
    nCells = shp[1]*shp[2]
    idxShuffle = sample(collect(1:1:nCells),nCells,replace=false)
    for i in idxShuffle
        cellIndex = ca.pa_cart_index[i]
        if ca.pa[i] === 1.0
            # Determine the neighbourhood weighted survival probability
            neighbourWeight = neighbourHoodWeight(ca.pa,env.suitabilityActive,cellIndex,shp,ca.neighbourWeightMatrix,ca.r)
            survWeight = neighbourWeight * ca.neighbourSurvivalWeight
            # Logistic function to scale survival proba
            sf = 1/(1+MathConstants.e^(-10*(env.suitabilityActive[cellIndex]-(1-survWeight))))
            survivalProbability = env.suitabilityActive[cellIndex] + (1-env.suitabilityActive[cellIndex])*sf
            # Determine survival
            survived = sample(rng,[0,1],ProbabilityWeights([1.0-survivalProbability,survivalProbability]))
            if survived === 0
                ca.pa[cellIndex] = 0
            end
        end
    end
end

end # Module
