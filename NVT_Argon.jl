using LinearAlgebra, Plots

mutable struct p
    m
    x
    c
    function p(m, x, c=0) 
        return new(m, x, c)
    end
end


function LJ(p1, p2, r_cut)
    if ((p1.x - p2.x)'*(p1.x - p2.x))[1] > r_cut^2
        return 0
    end
    r = norm(p1.x - p2.x)
    return 4*r^(-6)*(r^(-6) - 1)
end

function energy(p, particles)
    result = 0
    for i in particles
        if i != p
            result += LJ(p, i, r_cut)
        end
    end
    return result - (8/3)*π*((1/3)*(r_cut^(-9) - r_cut^(-3)))
end

function move(particles,dx)
    o = rand()*length(particles) |> ceil |> Int
    e_old = energy(particles[o], particles)

    p_new = p(m, particles[o].x + (rand(3,1)-0.5*ones(3,1))*dx, 1)
    e_new = energy(p_new, vcat(particles[1:o-1], particles[o+1:end]))
    kb = 1.38e-23
    t = exp(-kb^(-1)*(e_new - e_old))
    if rand() < t
        particles = replace(x -> x == particles[o] ? p_new : x, particles)
    end
    return particles
end

function run(particles, dx, iter, sample=100)
    result = zeros(1, div(iter, sample) |> floor |> Int)
    d = 0
    for i = 1:iter
        particles = move(particles,dx)
        if i % sample == 0
            d += 1
            for j = 1:length(particles)
                for k = 1:length(particles)
                    if j < k
                        result[d] = LJ(particles[j], particles[k], r_cut)
                    end
                end
            end
        end
    end
    return result
end


β, m = 1/(1.38e-23), 1
r_cut = 2.5

result = zeros(5,1000)
for i = 1:5
    particles = []
    for j = 1:10
        push!(particles, p(m, 10^(1/3)*rand(3,1)))
    end
    result[i,:] = run(particles, 0.005, 1e5)
end



plot(result'[100:500,1:4], title="N=10, iter=1e4")
savefig("./argon.pdf")

print(length(e))
σ = 3.405e-10
rho = (1680*σ^3)/.03994*6.023e23

