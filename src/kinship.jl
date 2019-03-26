###############################################################
# function to calculate kinship from genotype probability array
###############################################################

function calcKinship(geno::Matrix{Float64})

    # get dimensions
    sz = size(geno)

    # assign to variables for convenience
    nr = sz[1]
    nc = sz[2]

    # if empty then there is nothing to do
    if(nr==0)
        error("Nothing to do here.")
    else
        # make matrix to hold distances
        d = zeros(nr,nr)
    end

    # assign diagonals to ones
    for i=1:nr
        d[i,i] = 1.0
    end

    ncomplete = nc        
    # off-diagonal elements need to be calculated    
    if(nr>=2)
        for i=1:(nr-1)
            for j=(i+1):nr
                p1 = g[i,:]
                p2 = g[j,:]
                d[i,j] = d[j,i] = sum( p1 .* p2
                                       + (1 .- p1) .* (1 .- p2) ) / ncomplete
                
            end
        end
    end
    return d
end


function calcKinship(geno::Matrix{Union{Missing,Float64}})

    # get dimensions
    sz = size(geno)

    # assign to variables for convenience
    nr = sz[1]
    nc = sz[2]

    # if empty then there is nothing to do
    if(nr==0)
        error("Nothing to do here.")
    else
        # make matrix to hold distances
        d = zeros(nr,nr)
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
                iscomplete = .!( ismissing.(g[i,:]) .& ismissing.(g[j,:]) )
                ncomplete = sum(iscomplete)
                p1 = g[i,iscomplete]
                p2 = g[j,iscomplete]
                d[i,j] = d[j,i] = sum( p1 .* p2
                                       + (1-p1) .* (1-p2) ) / ncomplete
                
            end
        end
    end
    return d
end

