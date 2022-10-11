using LinearAlgebra
using PhysicalConstants: CODATA2018
using Unitful
using Distributions
using SparseArrays

using Unitful: eV, nm

abstract type QuantumSystem end

"""
    QuantumWire

Object holding all static information about the quantum wire system.

# Fields

| Field/Attribute  |   Type   |                                            Description                                                             |
|:----------------:|:---------|:-------------------------------------------------------------------------------------------------------------------| 
| Uix              | Matrix   | Transformation matrix from the uncoupled photon-molecule states to the polariton basis (eigenvector matrix of H).  |
| evals            | Vector   | Vector containing the eigenvalues of the system.                                                                   |
| phot_range       | Range    | Range of pure photonic states in the uncoupled basis.                                                              | 
| mol_range        | Range    | Range of pure molecular states in the uncoupled basis.                                                             | 
| mol_energies     | Vector   | Vector containing molecular excitation energies.                                                                   |
| mol_positions    | Vector   | Vector containing molecular positions along the wire.                                                              |
| phot_energies    | Vector   | Vector containing energies for the cavity/wire radiation modes.                                                    |
"""
struct QuantumWire{R,C} <: QuantumSystem
    Uix::Matrix{C}
    evals::Vector{R}
    phot_range::UnitRange{Int64}
    mol_range::UnitRange{Int64}
    mol_energies::Vector{R}
    mol_positions::Vector{R}
    phot_energies::Vector{R}
    phot_wavevectors::Vector{R}
end

"""
    QuantumWire

Constructor for the object `QuantumWire` which stores all relevant static information for the system.

Molecules are organized in one dimension and are assumed to have an uniform dipole moment. 
Arguments of this function are taken with units (Quantity objects) from the `Unitful` packaged. For example:

```julia
julia> using Unitful: eV
julia> x = 1eV
1 eV
julia> unit(x)
eV
```

# Arguments

> Note: All arguments are KEYWORD arguments; thus, the function must be called as `build_hamiltonian(ΩR = 0.3eV, ϵ = 3, ...)`

| Argument  |   Type   |                                       Description                                                 |
|:---------:|:---------|:--------------------------------------------------------------------------------------------------| 
| ΩR        | Quantity | Rabi frequency: ΩR = μ₀√(ħω₀ρ/2ϵ).                                                                |
| ϵ         | Number   | Index of refraction squared: ϵ = n^2.                                                             |
| Nc        | Int      | Number of positive cavity modes. The total number of cavity modes is Ntotal = 2Nc + 1.            |
| Nm        | Int      | Number of molecules in the wire.                                                                  |
| a         | Quantity | Average distance between molecules.                                                               |
| σa        | Number   | Standard deviation w.r.t to the regular positions. Units assumed to match `a`.                    |
| ωM        | Quantity | Average molecular excitation energy. Supports frequency, wavenumber as well as true energy values.|
| σM        | Number   | Standard deviation for molecular excitation energies. Units assumed to match `ωM`.                |
| Ly        | Quantity | Cavity/wire length along Y direction. (Mirror).                                                   |
| Lz        | Quantity | Cavity/wire length along Z direction. (Mirror).                                                   |
| nz        | Int      | Radiation quantum number associated with Z coordinate.                                            |
| ny        | Int      | Radiation quantum number associated with Y coordinate.                                            |
"""
function QuantumWire(;ΩR::Quantity, ϵ::Number, Nc::Int, Nm::Int, a::Quantity, σa::Number, ωM::Quantity, σM::Number, 
    Ly::Quantity, Lz::Quantity, nz::Int, ny::Int, printout=false)

    ### 1. HANDLE UNITS AND PHYSICAL CONSTANTS

    # Get physical constant in standard units
    h = CODATA2018.PlanckConstant
    c = CODATA2018.SpeedOfLightInVacuum

    # Include units for the sigmas (for conversion purpuses)
    σa = σa * unit(a)
    σM = σM * unit(ωM)

    # Log
    if printout
        _print_inputs(ΩR, Nm, a, σa, ωM, σM, Nc, Ly, Lz, ny, nz, ϵ)
    end

    # Check input units for Rabi splitting and molecular excitation
    # If they are given in angular frequency (ω) or wavenumber (ν) units, 
    # convert it to energy using E = ħω = hcν
    invM = dimension(u"1m^-1")
    invS = dimension(u"1s^-1")

    # Convert units given in wavenumber to proper energy
    if dimension(ΩR) == invM
        ΩR = h*c*ΩR

    # Convert frequency to energy
    elseif dimension(ΩR) == invS
        ΩR = h*ΩR/2π  # Note, we assume the frequency is angular
    end

    # Convert units given in wavenumber to proper energy
    if dimension(ωM) == invM
        ωM = h*c*ωM
        σM = h*c*σM

    # Convert frequency to energy
    elseif dimension(ωM) == invS
        ωM = h*ωM/2π  # Note, we assume the frequency is angular
        σM = h*σM/2π
    end

    # Strip units
    _ΩR = ustrip(u"eV", ΩR)
    _ωM = ustrip(u"eV", ωM)
    _σM = ustrip(u"eV", σM)
    _a  = ustrip(u"nm", a)
    _σa = ustrip(u"nm", σa)
    _Lz = ustrip(u"nm", Lz)
    _Ly = ustrip(u"nm", Ly)
    _ħ  = ustrip(u"eV*ps", h) / 2π
    _c  = ustrip(u"nm/ps", c)

    ### 2. ADD DISORDER

    # Create an array of energy values for molecules using a normal distribution 
    ωMvals = rand(Normal(_ωM, _σM), Nm)

    # Create an array with molecular positions uniformily spaces by `a` nm 
    # Plus a random deviation sampled from a normal distribution
    @assert Nm > 1

    avals = zeros(Nm)

    for i = 2:Nm # Note that the first molecule is fixed at x = 0
        avals[i] = (i-1)*_a + rand(Normal(0, _σa))
    end
    # Make sure positions are sorted
    if !issorted(avals)
        sort!(avals)
    end

    ### 3. PREPARE CAVITY MODES

    # Compute the wavenumber component associated with y and z directions.
    q₀ = √((nz*π/_Lz)^2 + (ny*π/_Ly)^2)

    # Length along x
    Lx = _a*Nm
    if printout
        output("\nLength along x (Nm * a) {:10.8f} nm", Lx)
    end

    # Prepare a vector with cavity energies and wavevectors. 
    Ncm = 2*Nc + 1
    wvec = zeros(Ncm)
    ωc = zeros(Ncm)
    wvec[1] = 0.0       # Minimal wavevector qx = 0
    ωc[1]   = _ħ*_c*q₀/√ϵ # Energy for qx = 0
    n = 2
    for mₓ = 1:Nc
        qₓ = 2π * mₓ / Lx
        ω = _ħ*_c*√((q₀^2 + qₓ^2)/ϵ)

        # Must include positive and negative components
        wvec[n] = qₓ
        wvec[n+1] = -qₓ

        # Energy is degenerate for -q and +q values
        ωc[n] = ω
        ωc[n+1] = ω
        n = n + 2
    end

    if printout
        output("\n Minimum Cavity energy: {:15.10f}", ωc[1])
        output(" Maximum Cavity energy: {:15.10f}", ωc[end])
        output("\n # Total System Dimensions: {}x{}", Nm+Ncm, Nm+Ncm)
        output("Memory required for Hamiltonian Matrix: {:5.3f} Gb | {:5.3f} Mb", 16e-9*(Nm+Ncm)^2, 16e-6*(Nm+Ncm)^2)
    end

    ### 4. CREATE HAMILTONIAN

    # Call unitless Hamiltonian converting input units to presumed ones
    H = build_hamiltonian(_ΩR, ωc, wvec, avals, ωMvals, printout=printout) 

    printout ? output("\n## Diagonalization") : nothing
    printout ? output("BLAS Config: {}", BLAS.get_config()) : nothing
    printout ? output("BLAS Number of Threads: {}", BLAS.get_num_threads()) : nothing
    printout ? output("\nDiagonalization started... ", ending="") : nothing
    t = @elapsed begin
        e, Uix = eigen(H)
    end
    printout ? output("Done in {:10.5f} seconds.", t) : nothing


    ### 5. RETURN SYSTEM OBJECT

    # We assume the photon basis come first, thus phot_range = 1:Np (Np being the number of photon states)
    # which in turn is Np = 2Nc + 1 (Nc is the number of positive momentum photonic states)
    # Similarly, the molecular range starts at Np+1 and ends at Np + Nm, where Nm is the number of molecules
    # That is 2Nc + 2 -> 2Nc + Nm + 1
    # This is the range of pure photonic states in the basis used to construct the Hamiltonian
    rphot = 1:(2Nc+1)
    # This is the range of pure molecular states in the basis used to construct the Hamiltonian
    rmol = (2Nc+2):(2Nc+Nm+1)

    # Construct and return the `QuantumWire` system object.
    return QuantumWire(Uix, e, rphot, rmol, ωMvals, avals, ωc, wvec)
end

struct SQuantumWire{R,C} <: QuantumSystem
    H::SparseMatrixCSC{C,Int64}
    phot_range::UnitRange{Int64}
    mol_range::UnitRange{Int64}
    mol_energies::Vector{R}
    mol_positions::Vector{R}
    phot_energies::Vector{R}
    phot_wavevectors::Vector{R}
end

function SQuantumWire(;ΩR::Quantity, ϵ::Number, Nc::Int, Nm::Int, a::Quantity, σa::Number, ωM::Quantity, σM::Number, 
    Ly::Quantity, Lz::Quantity, nz::Int, ny::Int, printout=false)

    ### 1. HANDLE UNITS AND PHYSICAL CONSTANTS

    # Get physical constant in standard units
    h = CODATA2018.PlanckConstant
    c = CODATA2018.SpeedOfLightInVacuum

    # Include units for the sigmas (for conversion purpuses)
    σa = σa * unit(a)
    σM = σM * unit(ωM)

    # Log
    if printout
        _print_inputs(ΩR, Nm, a, σa, ωM, σM, Nc, Ly, Lz, ny, nz, ϵ)
    end

    # Check input units for Rabi splitting and molecular excitation
    # If they are given in angular frequency (ω) or wavenumber (ν) units, 
    # convert it to energy using E = ħω = hcν
    invM = dimension(u"1m^-1")
    invS = dimension(u"1s^-1")

    # Convert units given in wavenumber to proper energy
    if dimension(ΩR) == invM
        ΩR = h*c*ΩR

    # Convert frequency to energy
    elseif dimension(ΩR) == invS
        ΩR = h*ΩR/2π  # Note, we assume the frequency is angular
    end

    # Convert units given in wavenumber to proper energy
    if dimension(ωM) == invM
        ωM = h*c*ωM
        σM = h*c*σM

    # Convert frequency to energy
    elseif dimension(ωM) == invS
        ωM = h*ωM/2π  # Note, we assume the frequency is angular
        σM = h*σM/2π
    end

    # Strip units
    _ΩR = ustrip(u"eV", ΩR)
    _ωM = ustrip(u"eV", ωM)
    _σM = ustrip(u"eV", σM)
    _a  = ustrip(u"nm", a)
    _σa = ustrip(u"nm", σa)
    _Lz = ustrip(u"nm", Lz)
    _Ly = ustrip(u"nm", Ly)
    _ħ  = ustrip(u"eV*ps", h) / 2π
    _c  = ustrip(u"nm/ps", c)

    ### 2. ADD DISORDER

    # Create an array of energy values for molecules using a normal distribution 
    ωMvals = rand(Normal(_ωM, _σM), Nm)

    # Create an array with molecular positions uniformily spaces by `a` nm 
    # Plus a random deviation sampled from a normal distribution
    @assert Nm > 1

    avals = zeros(Nm)

    for i = 2:Nm # Note that the first molecule is fixed at x = 0
        avals[i] = (i-1)*_a + rand(Normal(0, _σa))
    end
    # Make sure positions are sorted
    if !issorted(avals)
        sort!(avals)
    end

    ### 3. PREPARE CAVITY MODES

    # Compute the wavenumber component associated with y and z directions.
    q₀ = √((nz*π/_Lz)^2 + (ny*π/_Ly)^2)

    # Length along x
    Lx = _a*Nm
    if printout
        output("\nLength along x (Nm * a) {:10.8f} nm", Lx)
    end

    # Prepare a vector with cavity energies and wavevectors. 
    Ncm = 2*Nc + 1
    wvec = zeros(Ncm)
    ωc = zeros(Ncm)
    wvec[1] = 0.0       # Minimal wavevector qx = 0
    ωc[1]   = _ħ*_c*q₀/√ϵ # Energy for qx = 0
    n = 2
    for mₓ = 1:Nc
        qₓ = 2π * mₓ / Lx
        ω = _ħ*_c*√((q₀^2 + qₓ^2)/ϵ)

        # Must include positive and negative components
        wvec[n] = qₓ
        wvec[n+1] = -qₓ

        # Energy is degenerate for -q and +q values
        ωc[n] = ω
        ωc[n+1] = ω
        n = n + 2
    end

    if printout
        output("\n Minimum Cavity energy: {:15.10f}", ωc[1])
        output(" Maximum Cavity energy: {:15.10f}", ωc[end])
        output("\n # Total System Dimensions: {}x{}", Nm+Ncm, Nm+Ncm)
        output("\n # Non-Zero Elements:       {}", Nm+Ncm + 2*Nm*Ncm)
        output("Memory required for Hamiltonian Matrix: {:5.3f} Gb | {:5.3f} Mb", 16e-9*(Nm+Ncm)^2, 16e-6*(Nm+Ncm)^2)
    end

    ### 4. CREATE HAMILTONIAN

    # Call unitless Hamiltonian converting input units to presumed ones
    H = build_sparse_hamiltonian(_ΩR, ωc, wvec, avals, ωMvals, printout=printout) 


    ### 5. RETURN SYSTEM OBJECT

    # We assume the photon basis come first, thus phot_range = 1:Np (Np being the number of photon states)
    # which in turn is Np = 2Nc + 1 (Nc is the number of positive momentum photonic states)
    # Similarly, the molecular range starts at Np+1 and ends at Np + Nm, where Nm is the number of molecules
    # That is 2Nc + 2 -> 2Nc + Nm + 1
    # This is the range of pure photonic states in the basis used to construct the Hamiltonian
    rphot = 1:(2Nc+1)
    # This is the range of pure molecular states in the basis used to construct the Hamiltonian
    rmol = (2Nc+2):(2Nc+Nm+1)

    # Construct and return the `QuantumWire` system object.
    return SQuantumWire(H, rphot, rmol, ωMvals, avals, ωc, wvec)
end