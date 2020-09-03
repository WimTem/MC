include("./modules/Particle.jl")
include("./modules/MC.jl")
using .Particle, .MC
#Argon
σ = 3.405e-10   #[m[]]
M = 0.03994    #[kgMol-1]
#T = 1          <-> T = 119.8K
# rho = 1.0     <-> rho = 1680 [kgm-3]
# dt = 0.005    <-> dt = 1.09e-14 [s]
# P = 1         <-> P = 41.9 [MPa]
m = M/(6.022e23)


β = 1/(1.38e-23)
particles = []
for i = 1:8
    for j = 1:8
            push!(particles, Particle.p(m, σ.*[i, j]))
    end
end

result = MC.run(particles, 10, 2.5, 0.005)


using Plots

plot(result)