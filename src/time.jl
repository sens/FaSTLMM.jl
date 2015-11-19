using DataFrames

function runtime(funlist,n)
    # number of functions to be compared
    nf = length(funlist)

    elapsedTimes = Array(Float64,(n,nf))
    medianTimes = Array(Float64,nf)
    idx = 1

    for i=1:n
        for j=1:nf
            t0 =time()
            eval(parse(funlist[j]))
            t1 = time()
            elapsedTimes[i,j] = t1-t0
        end
    end

    for j in 1:nf
        medianTimes[j] = median(elapsedTimes[:,j])
    end

    return elapsedTimes, medianTimes
end

