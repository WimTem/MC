module Particle
    mutable struct p
        m
        x
        function p(m, x) 
            return new(m, x)
        end
    end
end