
##############################################################
# Routines for reading in genotypes and phenotypes
##############################################################

# For now, we are reading in data from Pjotr Prin's format.  Other
# formats will be included later, as needed.

# using DataArrays

using DelimitedFiles

function readPheno(file::AbstractString,nSkip::Int64,
                   nPheno::Int64,nInd::Int64)

    # allocate space for phenotypes
    pheno = Array{Union{T, Missing}}(Float64,nInd,nPheno)
    
    f = open(file,"r")
        
    for i=1:nSkip
        aLine = readline(f)
    end

    for i=1:nInd
        aline = readline(f)
        words = split(aline)
        for j=1:nPheno
            pheno[i,j] = str2num(words[j+1])
        end
    end
    close(f)
        
    return pheno
end


#####################################################################        

# function assumes that the first line is marker names, and the first
# column is ids; both are strings
#
# the rest of the file is assumed to be probablities, i.e. numeric

function readGenoProb(file::AbstractString;
                      dlm::AbstractChar=',',
                      getmarkernames::Bool=true,
                      getids::Bool=true)

    if(getmarkernames)
        # read file
        d = readdlm(file,dlm,header=true)
        # markernames
        markernames = d[2][2:end]
        # genotype probabilities
        gd = d[1]
    else
        # read file
        d = readdlm(file,dlm=dlm,header=true)
        # genotype probabilities
        gd = d
    end

    # get ids        
    if(getids)
        ids = convert(Vector{String},gd[:,1])
        gp = convert(Matrix{Float64},gd[:,2:end])
        markernames = markernames[1:end-1]
    else
        gp = convert(Matrix{Float64},gd)
    end

    # returns only genotype probabilities; others discarded for now    
    return gp
end


    
#####################################################################        
function readGeno(file::AbstractString,nSkip::Int64,
                  nMarkers::Int64,nInd::Int64,format::AbstractString="HAB")

    # allocate space for marker names and genotypes
    mNames = Array{String,1}(undef, nMarkers)
    geno = Array{Float64,2}(undef, nInd,nMarkers)
    
    if(format!="HAB")
        error("Cannot read this type of format.")
    end

    f = open(file,"r")
        
    for i=1:nSkip
        aLine = readline(f)
    end

    for i=1:nMarkers
        aline = readline(f)
        words = split(aline)
        mNames[i] = words[1]
        geno[:,i] = word2array(words[2],nInd)
    end
    close(f)
        
    return mNames, geno
end

#####################################################################    
function word2array(word::SubString{String},wordLen::Int64)
    g = Array{Int64,1}(undef, wordLen)
    for i=1:wordLen
        g[i] = f2code(word[i])
    end
    return g
end

#####################################################################    
function f2code(x::Char)
    return if(x=='A')
             0
           elseif(x=='H')
             1
           elseif(x=='B')
             2
           else
             NA
           end
end        

#####################################################################    
function str2num(x::SubString{String})
    n = tryparse(Float64,x)
    return isnull(n) ? NA : get(n)
end


    
