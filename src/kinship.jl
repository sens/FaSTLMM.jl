###############################################################
# function to calculate kinship from genotype probability array
###############################################################

using Distances

function calcKinship(geno::Array{Float64,2})

    # get dimensions
    sz = size(geno)

    # assign to variables for convenience
    nr = sz[1]
    nc = sz[2]

    # make matrix to hold distances
    d = Array{Float64,2}(undef,nr,nr)

    # if empty then there is nothing to do
    if(nr==0)
        error("Nothing to do here.")
    end

    # assign diagonals to ones
    for i=1:nr
        d[i,i] = 1.0
    end

    iscomplete = Array{Bool,1}(undef,nc)
    ncomplete::Int64 = 0    
    # off-diagonal elements need to be calculated    
    if(nr>=2)
        for i=1:(nr-1)
            for j=(i+1):nr
                iscomplete = !. ( ismissing.(g[i,:]) &. ismissing.(g[j,:]) )
                ncomplete = sum(iscomplete)
                d[i,j] = d[j,i] = sum(g[i,iscomplete]*g[j,iscomplete] +
                                      (1-g[i,iscomplete])*(1-g[j,iscomplete]))/(ncomplete)
                
            end
        end
    end
    return d
end

