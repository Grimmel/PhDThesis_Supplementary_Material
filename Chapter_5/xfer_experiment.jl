include("ca.jl")
using .CellularAutomata
using DelimitedFiles
using DataFrames
using CSV
using GaussianRandomFields

cal_landscapes = [787,663,686,218,168,753,36,212,576,80,777,328,547,682,90,456,688,224,775,321]
val_landscapes =[32,687,769,103,812]
function idwDistWeightMatrix(r,scale)
    w = 2*r+1
    m = ones(w,w)
    mIdx = CartesianIndices(m)
    for idx in mIdx
        #print(idx[1])
        m[idx] = sqrt(((idx[1]-(r+1))^2)+((idx[2]-(r+1))^2))
    end
    m = (m.^-1)/scale
    m[r+1,r+1] = 0.0
    return(m)
end

function rescale(x,newMin,newMax,oldMin,oldMax)
    return (newMax-newMin).*((x.-oldMin)./(oldMax-oldMin)).+newMin
end
function simulateNoNB(ls)
    suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability"*string(ls)*".asc",skipstart=6)
    cov = CovarianceFunction(2, Matern(1/8,0.5, σ=0.06,p=2))
    pts = range(0, stop=1, step=1/399)
    grf = GaussianRandomField(0.0,cov, GaussianRandomFields.CirculantEmbedding(), pts, pts,minpadding = 100)
    params = [["do1",1.0,0.2,2],["do2",5.0,0.4,5]]
    for p in params
        name = p[1]
        meanDisp = p[2]
        dispProb = p[3]
        maxDisp = p[4]
        for i in 1:10
            simName = name*"_ls"*string(ls)*"_rep"*string(i)
            pa = ones(400,400)
            paIdx = CartesianIndices(pa)
            ca = OccurenceCellularAutomata(pa,paIdx,deepcopy(suit),meanDisp,dispProb,maxDisp)
            for i in 1:100
                ca.suitability = deepcopy(suit) .+ sample(grf)
                for i in 1:400
                    for j in 1:400
                        if ca.suitability[i,j] <0
                            ca.suitability[i,j] = 0
                        elseif ca.suitability[i,j] >1
                            ca.suitability[i,j] = 1
                        end
                    end
                end
                coloniseSuitWeight(ca)
                extinction(ca)
            end
            open("D:/PHDExperimentOutputs/Transferability/occurrence/"*simName*"_es.asc","w") do io
                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                writedlm(io,pa)
            end
        end
    end
end
#params = [["strong1",0.65,1.0,0.3,2,1.2,1.2,3],["strong2",0.6,5.0,0.8,5,1.2,1.2,3]]
# params = [["strong4",0.55,3.0,0.8,5,1.2,1.2,3],["strong5",0.7,3.0,0.7,5,1.2,1.2,3]]

function simulateNB(ls)
    #params = [["mod1",0.9,1.0,0.45,2,0.5,0.5,1],["mod2",0.85,5.0,0.4,5,0.5,0.5,1],["strong1",0.65,1.0,0.4,2,1.2,1.2,3],["strong2",0.6,5.0,0.4,5,1.2,1.2,3]]
    params = [["vstrong1",0.5,1.0,0.1,2,2.5,1.5,3],["vstrong2",0.5,5.0,0.1,5,2.5,2.5,3]]
    cov = CovarianceFunction(2, Matern(1/8,0.5, σ=0.06,p=2))
    pts = range(0, stop=1, step=1/399)
    grf = GaussianRandomField(0.0,cov, GaussianRandomFields.CirculantEmbedding(), pts, pts,minpadding = 100)
    for p in params
        name = p[1]
        scale = p[2]
        meanDisp = p[3]
        dispProb = p[4]
        maxDisp = p[5]
        nsw = p[6]
        ndw = p[7]
        r = p[8]
        weight = idwDistWeightMatrix(r,10)
        for i in 1:10
            simName = name*"_ls"*string(ls)*"_rep"*string(i)
            pa = ones(400,400)
            paIdx = CartesianIndices(pa)
            suit = readdlm("D:/PHDExperimentOutputs/Transferability/landscapes/suitability/suitability"*string(ls)*".asc",skipstart=6)
            suit = rescale(suit,0,scale,0,1)
            ca = OccurrenceCellularAutomataNB(pa,paIdx,deepcopy(suit),meanDisp,dispProb,maxDisp,nsw,ndw,weight,r)
            for i in 1:200
                ca.suitability = deepcopy(suit) .+ sample(grf)
                for i in 1:400
                    for j in 1:400
                        if ca.suitability[i,j] <0
                            ca.suitability[i,j] = 0
                        elseif ca.suitability[i,j] >1
                            ca.suitability[i,j] = 1
                        end
                    end
                end
                coloniseNeighbourWeight(ca)
                extinctionNeighbourWeight(ca)
            end
            open("D:/PHDExperimentOutputs/Transferability/occurrence/"*simName*"_es.asc","w") do io
                write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                writedlm(io,pa)
            end
        end
    end
end

function runSims(landscapeList)
    for ls in landscapeList
        simulateNB(ls)
    end
end

#runSims(val_landscapes)
runSims(val_landscapes)
