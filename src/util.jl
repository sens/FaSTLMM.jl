###################################
# utility functions
###################################

using Statistics
# using DataStructures

# centers each column
function colCenter!(A::Matrix{Float64})

    (n,m) = size(A)

    # get mean of each column; convert to vector    
    colMeans = mean(A,dims=1) |> vec
        
    for i=1:n
        for j=1:m
            A[i,j] = A[i,j] - colMeans[j]
        end
    end
end

# centers each row
function rowCenter!(A::Matrix{Float64})

    (n,m) = size(A)

    # get mean of each column; convert to vector    
    rowMeans = mean(A,dims=2) |> vec
        
    for i=1:m
        for j=1:n
            A[j,i] = A[j,i] - rowMeans[j]
        end
    end
end

function colDivide!(A::Matrix{Float64},x::Vector{Float64})

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

function colStandardize!(A::Matrix{Float64})

    colCenter!(A)
    s = std(A,dims=1) |> vec
    colDivide!(A,s)
    
end
    
function rowDivide!(A::Matrix{Float64},x::Vector{Float64})

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

function compareValues(x_true::Array{Float64,1}, x::Array{Float64,1}, tolerance::Float64, threshold::Float64)
    if size(x_true) != size(x)
        error("Dimention Mismatch! Must compare two arrays of same length!")
    end

    passes = falses(size(x_true))
    t_passes = falses(0)
    for i in 1:size(x_true)[1]
        e = abs(x[i]-x_true[i])
        if e <= tolerance
            passes[i] = true
        end
        
        if x[i] >= threshold
            if e <= tolerance 
                push!(t_passes, true)
            else
                push!(t_passes,false)
            end
        end
        
    end
    # pass_rate = sum(passes) / size(x_true)[1]
    # pass_rate = sum(t_passes) / size(t_passes)[1]
    
    return (sum(passes), size(t_passes)[1], sum(t_passes))

end



