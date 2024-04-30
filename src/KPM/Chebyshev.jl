using LinearAlgebra
using KrylovKit
using PolaritonicSystems: SymBlockArrowHead

"""
    cheb_Tn(H, n)

Given an input H returns the associated Chebyshev polynomial of order n T_n.

H can be a number, a matrix or any object for which the operations * and - with itself are defined.
"""
function cheb_Tn(H, n)

    if n == 0
        return one(H)
    elseif n == 1
        return H
    end

    Hᵢ = one(H)
    Hᵢ₊₁ = H
    for i in 2:n
        Hᵢ₋₁ = Hᵢ 
        Hᵢ = Hᵢ₊₁
        Hᵢ₊₁ = 2*H*Hᵢ - Hᵢ₋₁
    end

    return Hᵢ₊₁
end

"""
    eval_cheb(x, n)

Returns a vector of Tn(x) up to order n, where Tn(x) are the Chebyshev polynomials of order n evaluated at x.
"""
function eval_cheb(x, n)

    if n == 0
        return [1]
    end

    T = zeros(n+1)
    T[1] = 1
    T[2] = x
    for i in 3:(n+1)
        T[i] = 2*x*T[i-1] - T[i-2]
    end

    return T
end

"""
    trace_Tn(H, n)

Return an array of length n where each entry is the trace of Tn(H), where Tn is the n--th Chebyshev polynomial.
Note that, if H is a number, this function is equivalent to `eval_cheb`.

This function can be used to evaluate the Chebyshev moments of a density of states exacly. Hence, it is used here
for validation purpuses. 
"""
function trace_Tn(H, n)

    out = zeros(n+1)
    Hᵢ = one(H)
    Hᵢ₊₁ = H
    for i = 0:n
        if i == 0
            # For n = 0 Tn = I, thus (1/D)*Tr(I) = 1
            out[i+1] = size(H,1)
        elseif i == 1
            out[i+1] = real.(tr(H))
        else

            Hᵢ₋₁ = Hᵢ 
            Hᵢ = Hᵢ₊₁
            Hᵢ₊₁ = 2*H*Hᵢ - Hᵢ₋₁

            out[i+1] = real.(tr(Hᵢ₊₁))
        end
    end

    return out
end

function cheb_scale(H::AbstractArray{T}; ϵ=0.01) where T
    Emin = eigsolve(H, 1, :SR)[1][1]
    Emax = eigsolve(H, 1, :LR)[1][1]
    a = T.((Emax - Emin) / (2 - ϵ))
    b = T.((Emax + Emin) / 2)
    return (H - b*I) ./ a, a, b
end

function cheb_scale(H::SymBlockArrowHead{T}; ϵ=0.01) where T
    Emin = eigsolve(H, 1, :SR)[1][1]
    Emax = eigsolve(H, 1, :LR)[1][1]
    a = (Emax - Emin) / (2 - ϵ)
    b = (Emax + Emin) / 2
    newd = (H.d .- b) ./ a
    newX = H.X ./ a
    return SymBlockArrowHead(newd, newX, H.l, H.r1, H.r2), a, b
end

function cheb_undo_scale(v::Vector, a, b)
    return a .* v .+ b
end
