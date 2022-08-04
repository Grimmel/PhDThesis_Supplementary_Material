include("D:/git/SpatialVirtualSpecies/src/SpatialVirtualSpecies.jl")
using .SpatialVirtualSpecies
using DataFrames,DelimitedFiles
using PyPlot
using GaussianRandomFields
# SAC noise
cov = CovarianceFunction(2, Matern(1/8,0.5, σ=0.08, p=2))
pts = range(0, stop=1, step=1/399)
pts2 = range(0, stop=1, step=1/799)
grf = GaussianRandomField(0.0,cov, GaussianRandomFields.CirculantEmbedding(), pts, pts2,minpadding = 100)
# Dispersal
pos_params = SpatialVirtualSpecies.ExponentialPosSelector(10)
cParam = SpatialVirtualSpecies.ColoniseWeightedParameters(pos_params,0.15,3)
# Neighbourhood
nt = SpatialVirtualSpecies.MooreNeighbours(2)
survivalScaling = SpatialVirtualSpecies.LogisticScalingParameters()
coloniseScaling = SpatialVirtualSpecies.LogisticScalingParameters()
wtGen = SpatialVirtualSpecies.IDW_MooreNeighbourhoodWeight(2,1,true)
wtMatrix = SpatialVirtualSpecies.generateWeightMatrix(wtGen)
survWeight =  0.05
colWeight = 0.05
nbParams = SpatialVirtualSpecies.NeighbourhoodParameters(nt,survWeight,survWeight,survivalScaling,coloniseScaling,wtMatrix)
# Function to ensure suitability values remain between zero and one when adding SAC noise
function checkValueLimits(ca::Matrix{Float64})
    shp = size(ca)
    for i in 1:shp[1]        #ROWS
        for j in 1:shp[2]    #COLS
            if ca[i,j] <0
                ca[i,j] = 0
            elseif ca[i,j] >1
                ca[i,j] = 1
            end
        end
    end
end
# Simulation behaviour
function simulateNB(ca,cp,nb,iterations)
    for i in 1:iterations
        noise = sample(grf)
        ca.survivalControl_Active = ca.survivalControl_Passive .+ noise
        ca.survivalControl_Active = ca.survivalControl_Passive .+ noise
        checkValueLimits(ca.survivalControl_Active)
        checkValueLimits(ca.coloniseControl_Active)
        SpatialVirtualSpecies.colonise(ca,cp,nb)
        SpatialVirtualSpecies.extinction(ca,nb)
    end
end
# Run simulations for each landscapes
for ls in ["3","10","11","1","2","5","6","12","13"]
    header = SpatialVirtualSpecies.getHeader("F:/PhD/prediction_extent/data/ls/species2/LS"*ls*"_suit.asc",6)
    suitability = readdlm("F:/PhD/prediction_extent/data/ls/species2/LS"*ls*"_suit.asc",skipstart=6)
    for i in 1:500
        pa = SpatialVirtualSpecies.generateStateLayer(suitability,0.01,(0.4,1.0))
        paIdx = CartesianIndices(suitability)
        ca = SpatialVirtualSpecies.SpeciesCellularAutomata(pa,paIdx,suitability,suitability)
        simulateNB(ca,cParam,nbParams,300)
        open("F:/PhD/prediction_extent/sims/species2/LS"*ls*"_rep"*string(i)*".asc","w") do io
            write(io,header)
            writedlm(io,ca.pa)
        end
    end
end
