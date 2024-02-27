using LinearAlgebra
using Distributions
using FFTW

# Auxiliary function for testing purpuses
# Returns a eigenvalues and a random Hermitian matri
function _random_matrix(N; lmin=0, lmax=1)

    λ = lmin .+ lmax .* rand(N)

    x = rand(N,N)
    U = exp(-im * Hermitian(x))

    @assert U*U' ≈ I

    return λ, U * diagm(λ) * U'
end

function jackson_kernel(order::Int)
    N = order+1
    n = 0:(N-1)
    return ( (N .- n .+ 1) .* cos.(π .* n / (N+1)) + sin.(π .* n / (N+1)) ./ tan(π / (N+1)) ) / (N+1)
end

function get_Tn(H, n)

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

# Computes μn = (1/D)*Tr(Tn(H)) exacly for all μn up to n
# where D is the dimension of H , Tr is the trace operation, and Tn is the n-th Chebyshev polynomial
function exact_Tr_Tn(H, n)

    μs = zeros(n+1)
    Dinv = 1/size(H,1)
    Hᵢ = one(H)
    Hᵢ₊₁ = H
    for i = 0:n
        if i == 0
            # For n = 0 Tn = I, thus (1/D)*Tr(I) = 1
            μs[i+1] = 1
        elseif i == 1
            μs[i+1] = Dinv * real.(tr(H))
        else

            Hᵢ₋₁ = Hᵢ 
            Hᵢ = Hᵢ₊₁
            Hᵢ₊₁ = 2*H*Hᵢ - Hᵢ₋₁

            μs[i+1] = Dinv * real.(tr(Hᵢ₊₁))
        end
    end

    return μs
end

function brute_force_reconstruct(μs, xvals)

    out = similar(xvals)
    out .= 0.0

    for k in eachindex(xvals)
        xk = xvals[k]
        pf = 1/(π*sqrt(1-xk^2))
        Tn = eval_cheb(xk, length(μs)-1)
        out[k] = pf * (μs[1] + 2 .* dot(μs[2:end], Tn[2:end]))
    end

    return out
end

function get_DOS(H, xvals, order)

    μ = exact_Tr_Tn(H, order)
    gμ = μ .* jackson_kernel(order)

    return brute_force_reconstruct(gμ, xvals)
end

function get_DOS_μ(H, order, k)

    μs = zeros(ComplexF64, order+1)

    for _ in 1:k

        z = complex_unit(size(H,1))
        inner!(μs, H, z)
    end

    return real.(μs ./ (k*size(H,1)))
end

function get_local_DOS_μ(H, i::Int, order::Int)

    # Initialize α0 and α1 = H⋅α0
    α0 = similar(H, size(H, 1))
    α0 .= 0.0
    α0[i] = 1.0

    μs = similar(H, order+1)
    μs .= 0.0

    inner!(μs, H, α0)

    return  μs ./ size(H,1)
end

function get_local_DOS_μ(H, is::AbstractVector, order::Int)

    μs = zeros(eltype(H), order+1, length(is))
    μs .= 0.0
    α0 = similar(H, size(H, 1))

    for (n,i) in enumerate(is)
        # Initialize α0 and α1 = H⋅α0
        α0 .= 0.0
        α0[i] = 1.0

        inner!(@view(μs[:,n]), H, α0)
    end

    return  μs ./ size(H,1)
end

function brute_force_LDOS(H, order)

    μs = zeros(eltype(H), size(H,1), order+1)
    Dinv = 1/size(H,1)
    Hᵢ = one(H)
    Hᵢ₊₁ = H
    for i = 0:order
        if i == 0
            μs[:, i+1] .= Dinv
        elseif i == 1
            μs[:, i+1] .= Dinv * diag(H)
        else

            Hᵢ₋₁ = Hᵢ 
            Hᵢ = Hᵢ₊₁
            Hᵢ₊₁ = 2*H*Hᵢ - Hᵢ₋₁

            μs[:, i+1] .= Dinv * diag(Hᵢ₊₁)
        end
    end

    return μs
end

# Merge this functions with DOS?
function rademacher(n)
    return rand([-1.0, 1.0], n)
end

function complex_unit(n)
    dist = Uniform(0, 2π)
    return [exp(-im*θ) for θ = rand(dist, n)]
end

function inner!(μs, A, z0)

    order = length(μs) - 1

    z1 = A*z0
    ζ0 = z0⋅z0
    ζ1 = z0⋅z1

    μs[1] += ζ0
    μs[2] += ζ1
    μs[3] += 2 * (z1⋅z1) - ζ0

    zⱼ₋₂ = z0
    zⱼ₋₁ = z1
    for j = 2:ceil(Int, order/2)
        zⱼ = 2 * A*zⱼ₋₁ - zⱼ₋₂
        μs[2j] += 2 * (zⱼ₋₁⋅zⱼ) - ζ1 
        if 2j == order+1
            break
        end
        μs[2j+1] += 2 * (zⱼ⋅zⱼ) - ζ0 

        zⱼ₋₂ = zⱼ₋₁
        zⱼ₋₁ = zⱼ
    end

    return μs
end

function get_cheb_xk(M)
    return [cos(π *(k+0.5)/M) for k = 0:(M-1)]
end

# M is the number of points where the function will be evaluated over
function reconstruct(μs, M)

    γk = slow_dct(μs, M)
    xk = get_cheb_xk(M)

    for i in eachindex(γk)
        γk[i] = γk[i] / (π*sqrt(1-xk[i]^2))
    end

    return γk
end

function slow_dct(μs, M=length(μs))

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

function fast_dct(λs)

    # Get DCT type-II and undo normalization
    out = dct(λs) .* sqrt(length(λs)/2)
    out[1] = out[1] * sqrt(2) 

    return out
end

function integral_over_μ(g, μs, M)

    γk = slow_dct(μs, M)
    xk = get_cheb_xk(M)
    gk = [g(x) for x in xk]

    return (1/M) * (γk ⋅ gk)
end

function exact_Z(H, β)
    tr(exp(-β.*H))
end

function approx_Z(H, β, o)
    μs = get_DOS_μ(H, o, 20)
    mus = μs .* jackson_kernel(o)
    # adding plus one here because of the negative energies...
    return integral_over_μ(x->exp(-β*(x+1)), mus, 2*o)
end

