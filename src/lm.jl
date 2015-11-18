#############################################################
# Functions to perform linear regression with one response
#############################################################
#
# Notation:
# y = response
# X = design matrix
#
# Functions should work for matrix-valued y

# Uses elementary operators
function lm0(y,X)
    xx = X'X
    xy = X'y
    b = xx\xy
    return b
end

# Uses specialized function for multiplying transposed matrix
function lm1(y,X)
    xx = At_mul_B(X, X) 
    xy = At_mul_B(X, y) 
    b = xx\xy
    return b
end

# Uses QR decomposition, multiplication function
function lm2(y,X)
    (q,r) = qr(X)
    y = At_mul_B(q,y)
    xx = At_mul_B(r,r)
    xy = At_mul_B(r,y) 
    b = xx\xy
    return b
end

# Uses QR decomposition and uses R matrix directly
function lm3(y,X)
    (q,r) = qr(X)
    y = At_mul_B(q,y)
    b = r\y
    return b
end
