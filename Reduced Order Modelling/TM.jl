using Pkg
Pkg.add(["MAT", "LinearAlgebra", "FastGaussQuadrature"])
using MAT
using LinearAlgebra
using FastGaussQuadrature

function HerrorFreq(LTI1, LTI2, spectralinterval)
    A1, B1, C1 = LTI1
    A2, B2, C2 = LTI2

    nodes, weights = gausslegendre(1000)
    a, b = spectralinterval
        L2norm = 0.0

        Herr(s) = first(C1*((s*I-A1)\B1)) - first(C2*((s*I-A2)\B2))

    for (i, s) in enumerate(nodes)
        L2norm += weights[i] * Herr(((b - a) * s + b + a) / 2)^2
    end
    L2norm *= (b - a) / 2
    absHerr(s) = abs(Herr(s))
    return L2norm, absHerr
end

function truncatedmeasure(A, B, C, n; p = ())
    #TODO: Here we need to implement truncated measure (use algo1 or algo 2) 
    V = zeros(1000,n)
    for i in 1:n
        V[i,i] = 1
    end
    D,T = eigen(A, sortby = - )

    Vhat = T * V
    invT = inv(T)
    What = transpose(invT) * V

    Ahat = transpose(What) * A * Vhat
    Bhat = transpose(What)*B
    Chat = C*Vhat

    return (Ahat, Bhat, Chat)
end

let
    vars = matread("C:/Users/lukev/Projects/NSF_REU_2025/SMU_REU_2025/codes/fullLTI.mat")
    A = vars["A"]
    B = vars["B"]
    C = vars["C"]

    Â, B̂, Ĉ = truncatedmeasure(A, B, C, 50)
    println("size(A) = ", size(A))
    println("size(Â) = ", size(Â))
    lerr, _ = HerrorFreq((A, B, C), (Â, B̂, Ĉ), (1, 20.0))
    println("||H-Ĥ||ₕ = ", lerr)
end
