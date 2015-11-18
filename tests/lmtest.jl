# include the function
include("../code/lm.jl")

N = 1000
P = 1

o = ones(Float64,(N,1));
x = randn( RandomDevice(), (N,P) ) * 10 + 50 ;
e = randn( RandomDevice(), (N,1) ) *2.5 ;
b = [30, 0.5;]
X = [o x];
y = X*b + e;

@time lm0(y,X);
@time lm1(y,X);
@time lm2(y,X);
@time lm3(y,X);


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

