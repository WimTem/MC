using LinearAlgebra, Plots

T, N, steps = 10, 1000, 100000



function banana(T, N, steps)
    n = ones(N, 3)
    viz = ones(steps |> Int)
    ξ = 8.314e-2*T
    for i = 1:steps
        E = 3*N*ξ^2/2
        j = rand()*N |> ceil |> Int
        k = rand()*3 |> ceil |> Int
        if rand()<0.5
            dn = 1
            dE = (2*n[j,k]+1)*ξ^2/2
        else
            dn = -1
            dE = (-2*n[j,k]+1)*ξ^2/2
        end

        if n[j,k] > 1 || dn == 1
            if rand()<exp(-dE/T)
                n[j, k] += dn
                E += dE
            end
        end
        viz[i] = E
    end
    return viz, n
end

viz, n = banana(T, N, steps)

plot(viz[1:1000:end], label="Energy")


