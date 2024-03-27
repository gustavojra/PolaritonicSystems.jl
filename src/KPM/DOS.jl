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

#function get_full_LDOS(H, order; jackson=true)

    #Nsites = size(H,1)

    #out = zeros(order+1, Nsites)
    #for i in 1:Nsites
        #out[:, i] = get_LDOS(H, i, order, jackson=jackson)
    #end

    #return out
#end

function get_DOS_μ(H, order, k)

    μs = zeros(ComplexF64, order+1)

    for _ in 1:k

        z = complex_unit(size(H,1))
        inner!(μs, H, z)
    end

    return real.(μs ./ (k*size(H,1)))
end

function get_DOS_μ(H::SparseMatrixCSC, order, k)

    μs = [zeros(ComplexF64, order+1) for i in 1:Threads.nthreads()]

    Threads.@threads for _ in 1:k

        z = complex_unit(size(H,1))
        inner!(μs[Threads.threadid()], H, z)
    end

    return real.(sum(μs) ./ (k*size(H,1)))
end

function get_LDOS_μ(H, i::Int, order::Int)

    # Initialize α0 and α1 = H⋅α0
    α0 = sparsevec([i], [1.0], size(H,1))

    μs = zeros(eltype(H), order+1)
    μs .= 0.0

    inner!(μs, H, α0)

    return  μs ./ size(H,1)
end

# SLOW AF
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