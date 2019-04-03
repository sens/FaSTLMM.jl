###################################
# utility functions
###################################

function colScale!(A::Matrix{Float64},x::Vector{Float64})

    (n,m) = size(A)
    if(length(x)!=m)
        error("Matrix and vector size do not match.")
    end

    for i=1:n
        for j=1:m
            A[i,j] = A[i,j]/x[j]
        end
    end
end

function rowScale!(A::Matrix{Float64},x::Vector{Float64})

    (n,m) = size(A)
    if(length(x)!=n)
        error("Matrix and vector size do not match.")
    end

    for i=1:m
        for j=1:n
            A[j,i] = A[j,i]/x[j]
        end
    end
end

"""
perform random shuffles of vector
the first column is the original vector if original=true
"""
function shuffleVector(rng::AbstractRNG,x::Vector{Float64},
                 nshuffle::Int64,original::Bool=true)
    xx = zeros(length(x),nshuffle+1)
    xx[:,1] = x
    for i=1:nshuffle
        xx[:,i+1] = shuffle(rng,x)
    end

    return xx
end
