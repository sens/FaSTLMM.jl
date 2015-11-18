# include the function
include("../code/lm.jl")


# testing the lm function

N = 10000
P = 1000

X = randn( RandomDevice(), (N,P) );
e = randn( RandomDevice(), (N,1) );
b = [(1:P);]

y = X*b + e;

@time lm0(y,X);
@time lm1(y,X);
@time lm2(y,X);
@time lm3(y,X);

