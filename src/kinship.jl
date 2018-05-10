############################################
# function to calculate kinship
############################################

using Distances

function calcKinship(geno::DataArray{Float64,2})

    # get dimensions
    sz = size(geno)

    # assign to variables for convenience
    nr = sz[1]
    nc = sz[2]

    # make matrix to hold distances
    d = Array(Float64,nr,nr)

    # if empty then there is nothing to do
    if(nr==0)
        error("Nothing to do here.")
    end

    # assign diagonals to zero
    for i=1:nr
        d[i,i] = 0
    end

    iscomplete = Array(Bool,nc)
    ncomplete::Int64 = 0    
    # off-diagonal elements need to be calculated    
    if(nr>=2)
        for i=1:(nr-1)
            for j=(i+1):nr
                iscomplete = !(geno.na[i,:] | geno.na[j,:])
                ncomplete = sum(iscomplete)
                d[i,j] = d[j,i] = cityblock(geno.data[i,iscomplete],
                                      geno.data[j,iscomplete]) / (2*ncomplete)
                
            end
        end
    end
    return d
end

