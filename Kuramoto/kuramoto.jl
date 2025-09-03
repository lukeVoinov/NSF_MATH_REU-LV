using Random
using DifferentialEquations
using Distributions
using Plots
using MAT

function kuramoto(du, u, p, t)
    n, k, w = p

    for i = 1:n
        s = 0.0
        for j = 1:n
            s += k[i, j] * sin(u[j] - u[i])
        end
        du[i] = w[i] + s / n
    end
end

FT = Float64

f(theta) = (cos(theta), sin(theta))
colorcond(x) = 0 < mod(x, 2π) < π ? :red : :green
n = 10
fs = 250 #Sampling frequency
T = 1.0
theta0 = collect(0.0:2π/n:2π)
tspan = (0.0, T)
p = (N = n, K = 0.1 * ones(FT, (n, n)), w = rand(Uniform(0, 10), n))
prob = ODEProblem(kuramoto, theta0, tspan, p)
sol = solve(prob, RK4(), saveat = collect(0:1/fs:T), progress = true)

anim = Animation()
for t in sol.t
    scatter(
        f.(sol(t)),
        legend = false,
        xlims = (-1.1, 1.1),
        ylims = (-1.1, 1.1),
        aspect_ratio = 1.0,
        color = colorcond.(sol(t)),
    )
    frame(anim)
end
gif(anim, "kuramoto/kuramoto.gif", fps = 60)


matopen("kuramoto/kuramoto.data", "w") do file
    data = stack(sol.u);
    write(file, "data", data);
end

