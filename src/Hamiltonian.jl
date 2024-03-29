import Base: *
using LinearAlgebra

struct SymBlockArrowHead{T}
    d::Vector{T}
    X::Matrix{T}
    l::Int
    r1::UnitRange{Int64}
    r2::UnitRange{Int64}
end

function *(A::SymBlockArrowHead, v::Vector)

    if length(v) != A.l
        throw(DimensionMismatch("dimension of the vector [$(length(v))] does not match the required by the block arrowhead matrix [$(A.l)]"))
    end

    out = A.d .* v
    out[A.r2] .+= A.X * v[A.r1]
    out[A.r1] .+= A.X' * v[A.r2]

    return out
end

"""
    build_hamiltonian

Construct a Hamiltonian matrix under the Coulumb gauge and RWA approximation. 
See Ribeiro 2022: https://www.nature.com/articles/s42004-022-00660-0

Units are assumed to be eV for energy and nm for length. 

# Arguments
| Argument  |   Type   |                                       Description                                                 |
|:---------:|:---------|:--------------------------------------------------------------------------------------------------| 
| ΩR        | Quantity | Rabi frequency: ΩR = μ₀√(ħω₀ρ/2ϵ). In eV                                                          |
| ωc        | Vector   | Vector containing energies for the cavity/wire radiation modes.                                   |
| wvec      | Vector   | Vector containing wave-vector along x for the cavity/wire radiation modes.                        |
| avals     | Vector   | Vector containing molecular positions along the wire.                                             |
| ωMvals    | Vector   | Vector containing molecular excitation energies.                                                  |
"""
function build_hamiltonian(ΩR::T, ωc::Vector, wvec::Vector, avals::Vector, ωMvals::Vector; printout=false) where T <: AbstractFloat

    @assert length(ωMvals) == length(avals)

    ### 1. INITIALIZE HAMILTONIAN

    # Total size of the system
    Nm = length(ωMvals)
    Nc = length(ωc)
    Nt = Nm + Nc
    
    if printout
        output("\n## Constructing Hamiltonian Matrix")
        output("Allocating memory... "; ending="")
    end
    H = zeros(complex(T), Nt, Nt)
    printout ? output("Done.") : nothing

    ### 2. POPULATE DIAGONAL
    printout ? output("Populating diagonal... ", ending="") : nothing

    # Cavity energies
    for i = 1:Nc
        H[i,i] = ωc[i]
    end

    # Molecule energies
    for i = 1:Nm
        H[i+Nc, i+Nc] = ωMvals[i]
    end
    printout ? output("Done.") : nothing

    ### 3. POPULATE OFF-DIAGONAL
    printout ? output("Populating off-diagonal... ", ending="") : nothing

    # Fill off-diagonal terms with i < j
    # No need to fill j < i as we will declare the matrix Hermitian
    for i = 1:Nc
        for j = 1:Nm
            pf = -im * 0.5 * ΩR * √(ωMvals[j]/(Nm*ωc[i])) # LM interaction Prefactor
            x = avals[j]                                  # Position of the j-th molecule
            H[i,j+Nc] = pf * exp(-im * wvec[i] *x)
        end
    end
    printout ? output("Done.") : nothing

    return Hermitian(H)
end

function build_bah_hamiltonian(ΩR::Float64, ωc::Vector{Float64}, wvec::Vector{Float64}, avals::Vector{Float64}, ωMvals::Vector{Float64}; printout=false)
    if printout
        output("\n## Constructing Hamiltonian Matrix")
        output("Allocating memory... "; ending="")
    end
    d = zeros(ComplexF64, length(ωMvals)+length(ωc)) 
    X = zeros(ComplexF64, length(ωMvals), length(ωc)) 
    printout ? output("Done.") : nothing
    return build_bah_hamiltonian!(d, X, ΩR, ωc, wvec, avals, ωMvals, printout=printout)
end

function build_bah_hamiltonian(ΩR::Float32, ωc::Vector{Float32}, wvec::Vector{Float32}, avals::Vector{Float32}, ωMvals::Vector{Float32}; printout=false)
    if printout
        output("\n## Constructing Hamiltonian Matrix")
        output("Allocating memory... "; ending="")
    end
    d = zeros(ComplexF32, length(ωMvals)+length(ωc)) 
    X = zeros(ComplexF32, length(ωMvals), length(ωc)) 
    printout ? output("Done.") : nothing
    return build_bah_hamiltonian!(d, X, ΩR, ωc, wvec, avals, ωMvals, printout=printout)
end

function build_bah_hamiltonian!(d::Vector{T2}, X::Matrix{T2}, ΩR::T, ωc::Vector{T}, wvec::Vector{T}, 
                                avals::Vector{T}, ωMvals::Vector{T}; printout=false) where {T2 <: Complex, T <: AbstractFloat}

    @assert length(ωMvals) == length(avals)

    ### 1. INITIALIZE HAMILTONIAN

    # Total size of the system
    Nm = length(ωMvals)
    Nc = length(ωc)

    ### 2. POPULATE DIAGONAL
    printout ? output("Populating diagonal... ", ending="") : nothing
    # Cavity energies
    for i = 1:Nc
        d[i] = ωc[i]
    end

    # Molecule energies
    for i = 1:Nm
        d[i+Nc] = ωMvals[i]
    end
    printout ? output("Done.") : nothing

    ### 3. POPULATE OFF-DIAGONAL
    printout ? output("Populating off-diagonal... ", ending="") : nothing

    # Fill off-diagonal terms with i < j
    # No need to fill j < i as we will declare the matrix Hermitian
    for i = 1:Nm
        for j = 1:Nc
            pf = im * 0.5 * ΩR * √(ωMvals[i]/(Nm*ωc[j])) # LM interaction Prefactor
            x = avals[i]                                  # Position of the j-th molecule
            X[i,j] = pf * exp(-im * wvec[j] *x)
        end
    end
    printout ? output("Done.") : nothing

    return SymBlockArrowHead(d, X, Nc+Nm, 1:Nc, (Nc+1):(Nc+Nm))
end

function build_sparse_hamiltonian(ΩR::T, ωc::Vector, wvec::Vector, avals::Vector, ωMvals::Vector; printout=false) where T <: AbstractFloat

    @assert length(ωMvals) == length(avals)

    ### 1. INITIALIZE HAMILTONIAN

    # Total size of the system
    Nm = length(ωMvals)
    Nc = length(ωc)
    Nt = Nm + Nc

    Nnonzero = Nt + 2*Nc*Nm
    
    if printout
        output("\n## Constructing Sparse Hamiltonian Matrix")
    end

    Is = Int[]
    Js = Int[]
    # Create an array of values with data type matching ΩR, except that it is complex
    V = complex(T)[]

    ### 2. POPULATE DIAGONAL
    # Cavity energies
    for i = 1:Nc
        push!(Is, i)
        push!(Js, i)
        push!(V, ωc[i])
    end

    # Molecule energies
    for i = 1:Nm
        push!(Is, i+Nc)
        push!(Js, i+Nc)
        push!(V, ωMvals[i])
    end

    ### 3. POPULATE OFF-DIAGONAL

    # Fill off-diagonal terms with i < j
    # No need to fill j < i as we will declare the matrix Hermitian
    for i = 1:Nc
        for j = 1:Nm
            pf = -im * 0.5 * ΩR * √(ωMvals[j]/(Nm*ωc[i])) # LM interaction Prefactor
            x = avals[j]                                  # Position of the j-th molecule
            val = pf * exp(-im * wvec[i] *x)
            if abs2(val) > 1e-18
                push!(Is, i)
                push!(Js,j+Nc)
                push!(V, val)
                push!(Is,j+Nc)
                push!(Js, i)
                push!(V, conj(val))
            end
        end
    end
    printout ? output("Done.") : nothing

    return sparse(Is, Js, V, Nt, Nt)
end

# Experimental algorithm
# Testing new approach to leakage
function effstates_hamiltonian(ΩR::Number, ωc::Vector, wvec::Vector, avals::Vector, ωMvals::Vector, γ, Γ, nγ; printout=false)

    @assert nγ > 1

    # Total size of the system
    Nm = length(ωMvals)
    Nc = length(ωc)
    Nt = Nc + Nm + nγ

    H = zeros(ComplexF64, Nt, Nt)

    r = 1:(Nm+Nc)
    H[r,r] .= build_hamiltonian(ΩR, ωc, wvec, avals, ωMvals; printout=printout)

    printout ? output("Computing effective external cavity states") : nothing

    # Coupling between cavity modes and sink mode
    s = Nc+Nm+1 # Sink mode index
    for i = 1:Nc
        H[i,s] = γ
        H[s,i] = γ'
    end

    # Coupling between trap modes
    for i = s:(Nt-1)
        H[i,i+1] = Γ
        H[i+1,i] = Γ'
    end

    return Hermitian(H)
end