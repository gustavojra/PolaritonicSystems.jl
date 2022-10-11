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
function build_hamiltonian(ΩR::Number, ωc::Vector, wvec::Vector, avals::Vector, ωMvals::Vector; printout=false)

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
    H = zeros(ComplexF64, Nt, Nt)
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

function build_sparse_hamiltonian(ΩR::Number, ωc::Vector, wvec::Vector, avals::Vector, ωMvals::Vector; printout=false)

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
    V = ComplexF64[]

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