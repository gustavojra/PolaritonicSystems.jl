function DOS(H, order, k; jackson=true)

    mus = get_DOS_μ(H, order, k)

    if jackson
        apply_jackson_kernel!(mus)
    end

    return reconstruct(mus)
end

function LDOS(H, i, order; jackson=true)

    mus = get_LDOS_μ(H, i, order)

    if jackson
        apply_jackson_kernel!(mus)
    end

    return reconstruct(mus)
end

function LDOS(H, i, order, k; jackson=true)

    mus = get_LDOS_μ(H, i, order, k)

    if jackson
        apply_jackson_kernel!(mus)
    end

    return reconstruct(mus)
end

function get_DOS_μ(H::AbstractArray{T,2}, order, k) where T

    μs = zeros(T, order+1)

    for _ in 1:k
        z = T.(complex_unit(size(H,1)))
        inner!(μs, H, z)
    end

    return real.(μs ./ (k*size(H,1)))
end

function get_DOS_μ(H::SparseMatrixCSC{T}, order, k) where T

    μs = [zeros(T, order+1) for i in 1:Threads.nthreads()]

    Threads.@threads for _ in 1:k

        z = T.(complex_unit(size(H,1)))
        inner!(μs[Threads.threadid()], H, z)
    end

    return real.(sum(μs) ./ (k*size(H,1)))
end

function get_LDOS_μ(H::AbstractArray{T}, i::Int, order::Int) where T

    # Initialize α0 and α1 = H⋅α0
    α0 = sparsevec([i], [one(T)], size(H,1))

    μs = zeros(T, order+1)
    μs .= 0.0

    inner!(μs, H, α0)

    return  real.(μs ./ size(H,1))
end

function get_LDOS_μ(H::SymBlockArrowHead{T}, i::Int, order::Int) where T

    # Initialize α0 and α1 = H⋅α0
    α0 = zeros(T, size(H,1))
    α0[i] = 1.0

    μs = zeros(T, order+1)
    μs .= 0.0

    inner!(μs, H, α0)

    return  real.(μs ./ size(H,1))
end

function get_LDOS_μ(H::AbstractArray{T}, is::AbstractVector, order::Int, k::Int) where T

    μs = zeros(T, order+1)

    for _ in 1:k
        z = T.(complex_unit(size(H,1), is))
        inner!(μs, H, z)
    end

    return real.(μs ./ (k*size(H,1)))
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