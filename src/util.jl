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

