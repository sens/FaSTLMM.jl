###########################################################
# genome scan function; no covariates, two genotype groups
###########################################################

function scan(y::Array{Float64,2},g::Array{Float64,2},
              K::Array{Float64,2},reml::Bool)

    # number of markers
    (n,m) = size(g)
    # make intercept
    intcpt = ones(n,1)
    # rotate data
    (y0,X0,lambda0) = rotateData(y,[intcpt g],K)
    # fit null lmm
    out0 = flmm(y0,X0[:,1],lambda0)
    # weights proportional to the variances
    wts = makeweights( out0[:sigma2],out0[:h2],lambda0 )
    # rescale by weights
    scale!(sqrt.(1./wts),y0)
    scale!(sqrt.(1./wts),X0)

    # perform genome scan
    rss0 = sum(y0.^2)
    lod = zeros(m)
    for( i = 1:m )
        out = ls(y0,X0[:,i+1],reml,true)
        lod[i] = out[:ell]
    end

    return lod
        
end
              
function makeweights( sigma2::Float64, h2::Float64,
                      lambda::Array{Float64,1} )
    return h2*lambda+(1.0-h2)
end