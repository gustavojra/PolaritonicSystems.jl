module KPM
using LinearAlgebra
using Distributions
using FFTW
using KrylovKit
using SparseArrays

"""
    rademacher(n)

Returns a vector of length n with a Rademacher distribution, i.e. random entries of -1 and 1 with equal probability.
"""
function rademacher(n; T=Float64)
    return rand([-one(T), one(T)], n)
end

"""
    complex_unit(n)

Returns a vector where each entry is exp(-im*θ) where θ is draw from a uniform distribution [0, 2π]
"""
function complex_unit(n; T=Float64)
    dist = Uniform(0, 2π)
    return [exp(-im*θ) for θ = rand(dist, n)]
end

"""
    inner!(out, A, z0)    

Compute the inner product out[n] += ⟨z0|Tn(A)|z0⟩, where n refers to the order of the Chebyshev polynomial Tn(A). The order of the
expansion is taken as the length of the output vectot (out) minus one. Results are ADDED to the vector out (no overwritting).
"""
function inner!(out, A, z0)

    order = length(out) - 1

    z1 = A*z0
    ζ0 = z0⋅z0
    ζ1 = z0⋅z1

    out[1] += ζ0
    out[2] += ζ1
    out[3] += 2 * (z1⋅z1) - ζ0

    zⱼ₋₂ = z0
    zⱼ₋₁ = z1
    for j = 2:ceil(Int, order/2)
        zⱼ = 2 * (A*zⱼ₋₁) - zⱼ₋₂
        out[2j] += 2 * (zⱼ₋₁⋅zⱼ) - ζ1 
        if 2j == order+1
            break
        end
        out[2j+1] += 2 * (zⱼ⋅zⱼ) - ζ0 

        zⱼ₋₂ = zⱼ₋₁
        zⱼ₋₁ = zⱼ
    end

    return out
end

"""
    cheb_nodes(M)

Compute M Chebyshev nodes. These values are used for the polynomial interpolation using Chebyshev polynomials.
The nodes are equally spaced points of a unit circle projected onto the interval [-1, 1]. This allows one to used
the fast cosine transformation (see `reconstruct`) to compiute functions using Chebyshev moments.
"""
function cheb_nodes(M)
    return [cos(π *(k+0.5)/M) for k = 0:(M-1)]
end

"""
    reconstruct(μs)

Given a vector of Chebyshev moments μs, reconstruct the function over chebyshev nodes using a fast cosine transformation.
"""
function reconstruct(μs)

    γk = FFTW.r2r(μs, FFTW.REDFT01)
    xk = cheb_nodes(length(μs))

    for i in eachindex(γk)
        γk[i] = γk[i] / (π*sqrt(1-xk[i]^2))
    end

    return γk
end

"""
    _slow_dct(μs, M=length(μs))

Performs a discrete cosine transformation using a naive double loop algorithm O(N^2). Used for validation purpuses.
"""
function _slow_dct(μs, M=length(μs))

    out = zeros(M)

    for k = 0:(M-1)
        for i = eachindex(μs)
            n = i - 1

            if n == 0
                out[k+1] += μs[1] 
                continue
            end

            # Type-III DCT
            out[k+1] += 2 * μs[i] * cos(π*(k+1/2)*n/M)
        end
    end

    return out
end

"""
    jackson_kernel(order)

Computes the Jackson kernel coefficients that can be multiplied by Chebyshev cofficients.
"""
function jackson_kernel(order::Int)
    N = order+1
    n = 0:(N-1)
    return ( (N .- n .+ 1) .* cos.(π .* n / (N+1)) + sin.(π .* n / (N+1)) ./ tan(π / (N+1)) ) / (N+1)
end

"""
    apply_jackson_kernel!(μs)

Modify the input Chebyshev moments using the Jackson kernel.
"""
function apply_jackson_kernel!(μs)
    N = length(μs)

    for n in 0:(N-1)
        μs[n+1] *= ((N - n + 1) * cos(π*n/(N+1)) + sin(π*n/(N+1)) / tan(π/(N+1)) ) / (N+1)
    end
end

"""
    cheb_integral(g, μs, M)

Computes the integral ∫ f(x)⋅g(x) dx over the integral [-1, 1] where f is approximated using the Chebyshev moments μs. This integral
is performed in a trapezoidal manners, where g will be computed over the Chebyshev nodes and f is reconstructed over the same points.

∫ f(x)⋅g(x) dx ≈ (1/M) * ∑ₖ f(xk) g(xk)  where xk is a Chebyshev node.
"""
function cheb_integral(g, μs)

    γk = FFTW.r2r(μs, FFTW.REDFT01)
    M = length(γk)
    xk = cheb_nodes(M)
    gk = [g(x) for x in xk]

    return (1/M) * (γk ⋅ gk)
end

include("Chebyshev.jl")
include("DOS.jl")

end #module