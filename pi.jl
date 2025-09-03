using Pkg
Pkg.add(["Random","Statistics","Polynomials"])

using Random
using Statistics
using Polynomials

function getpi(n)
    x = rand(n)
    y = rand(n)
    c = count(xi^2 + yi^2 ≤ 1 for (xi,yi) in zip(x,y))
    p = 4*c/n
    return p
end

let
    numsteps = 3
    m = 50
    evals = zeros(numsteps)
    Ns = zeros(numsteps)

    for i = 1:numsteps
        N = 10^(2*i)
        pies = [getpi(N) for _ in 1:m]
        rmse = sqrt(mean((pies.-π).^2))
        Ns[i] = N
        evals[i] = rmse
    end

    logN = log10.(Ns)
    logE = log10.(evals)
    order = coeffs(fit(logN, logE, 1))[2]
    println(order)
end
