#######################################################################
# Functions to calculate residual sum of squares from linear regression
#######################################################################
#
# Notation:
# y = response
# X = design matrix
#
# Functions should work for matrix-valued y

# Uses elementary operators
function rss0(y,X)
    xx = X'X
    xy = X'y
    b = (xx\xy)
    yhat = X*b
    rss = norm(y-yhat)^2
    return rss
end

# Uses specialized function for multiplying transposed matrix
function rss1(y,X)
    xx = At_mul_B(X, X)
    xy = At_mul_B(X, y)
    b = xx\xy
    yhat = X*b
    rss = norm(y-yhat)^2
    return rss
end

# Uses QR decomposition, multiplication function
function rss2(y,X)
    (q,r) = qr(X)
    xy = At_mul_B(q,y)
    b = r\xy
    yhat = X*b
    rss = norm(y-yhat)^2
    return rss
end

# Uses full QR decomposition and computes residuals in transformed space
function rss3(y,X)
    (q,r) = qr(X,thin=false)
    d = size(r,1)
    z = At_mul_B(q,y)
    rss = norm(z[(d+1):end])^2
    return rss
end

# Uses extended QR decomposition and uses R matrix directly
function rss4(y,X)
    (q,r) = qr([X y])
    d = size(r,1)
    rss = r[d,d]^2
    return rss
end

# Uses Cholesky decomposition with augmented matrix
function rss5(y,X)
    Z = [X y]
    zz = Z'Z
    c = chol(zz)
    d = size(zz,1)
    rss = c[d,d]^2
    return rss
end

#

function rss6(y,X)

end
