using LinearAlgebra, Plots

ρ = 0.6
T = 2
β = ((1.38e-23)*(T))^(-1)
rc = 2.5
dt = 5e-3
σ = 3.405e-10

function U(p1, p2)
    r = norm(p1 - p2)
    return 4*((r)^(-12) - (r)^(-6))
end

function calc_energy(particles)
    result = 0
    for i in Iterators.product(1:size(particles)[1], 1:size(particles)[1]) |> collect
        if i[1] < i[2]
            result += U(particles[i[1], 1:2], particles[i[2], 1:2])
        end
    end
    return result
end

function move(particles, dt)
    U_old = calc_energy(particles)
    #Select random particle
    o = rand()*size(particles)[1] |> ceil |> Int
    #Add random displacement
    p_new = particles[o, 1:2] + (rand(2,1) - 0.5*ones(2,1))*dt
    #Periodic BC
    for i = 1:2
        if p_new[i] < 0
            p_new[i] += 1
        elseif p_new[i] > 1
            p_new[i] -= 1
        end
    end
    #Replace particle
    particles_temp = vcat(particles[1:o-1,:], p_new', particles[o+1:end,:])
    U_new = calc_energy(particles_temp)

    #Metropolis Criterion
    if (U_new - U_old) < 0
        return particles_temp
    elseif rand() < exp(-0.5*(U_new - U_old))
        return particles_temp
    end
    return particles
end


function run(particles, dt, iter, sample=500)
    energy = []
    for i = 1:iter
        particles = move(particles, dt)
        #Sample energy 
        if (i % sample) == 0
            push!(energy, calc_energy(particles))
        end
        #Check for convergence
        if length(energy) > 11 && abs(sum(energy[end-11:end-1])/10  - energy[end]) < 1
            print("\n Convergence in ", i, " iterations. \n")
            return energy, particles
        end
    end
    return energy, particles
end

particles = rand(16,2)

energy, particles2 = run(particles, 0.005, 200000)


scatter(particles[:,1], particles[:,2], label="Before")
scatter!(particles2[:,1], particles2[:,2], label="After")
savefig("./images/basic_particles.pdf")

plot(energy[100:end])
savefig("./images/basic_energy.pdf")

particles3 = rand(16,2)*1e-3 + 0.5*ones(16,2)
energy2, particles4 = run(particles, 0.005, 200000)

scatter(particles3[:,1], particles3[:,2], label="Before", xlims=(0,1), ylims=(0,1))
scatter!(particles4[:,1], particles4[:,2], label="After")
savefig("./images/basic_particles2.pdf")

plot(energy2[100:end])
savefig("./images/basic_energy2.pdf")