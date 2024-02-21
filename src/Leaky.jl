
struct LeakyQuantumWire{R,C} <: QuantumSystem
    Uix::Matrix{C}
    evals::Vector{R}
    phot_range::UnitRange{Int64}
    mol_range::UnitRange{Int64}
    mol_energies::Vector{R}
    mol_positions::Vector{R}
    phot_energies::Vector{R}
    phot_wavevectors::Vector{R}
    γ::R
    Γ::R
end

function LeakyQuantumWire(;ΩR::Quantity, ϵ::Number, Nc::Int, Nm::Int, a::Quantity, σa::Number, ωM::Quantity, σM::Number, 
    Ly::Quantity, Lz::Quantity, nz::Int, ny::Int, γ::Number, Γ::Number, nγ::Int, printout=false, balanced=true)

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

    ### 2. GET MOLECULAR SYSTEM
    ωMvals, avals = get_wire_molecules(Nm, _ωM, _σM, _a, _σa)

    ### 3. PREPARE CAVITY MODES

    # Length along x
    Lx = _a*Nm
    if printout
        output("\nLength along x (Nm * a) {:10.8f} nm", Lx)
    end

    # Get photon modes
    wvec, ωc = get_wire_modes(Nc, Lx, ny, _Ly, nz, _Lz, ϵ,  c=_c, ħ=_ħ)

    LeakyQuantumWire(_ΩR, ωc, wvec, ωMvals, avals, γ, Γ, nγ, printout=printout)
end

function LeakyQuantumWire(ΩR, ωc, wvec, ωMvals, avals, γ, Γ, nγ; printout=false)

    Ncm = length(wvec)
    Nm = length(ωMvals)

    if printout
        output("\n Minimum Cavity energy: {:15.10f}", ωc[1])
        output(" Maximum Cavity energy: {:15.10f}", ωc[end])
        output("\n # Total System Dimensions: {}x{}", Nm+Ncm, Nm+Ncm)
        output("Memory required for Hamiltonian Matrix: {:5.3f} Gb | {:5.3f} Mb", 16e-9*(Nm+Ncm)^2, 16e-6*(Nm+Ncm)^2)
    end

    ### CREATE HAMILTONIAN

    # Call unitless Hamiltonian converting input units to presumed ones
    H = effstates_hamiltonian(ΩR, ωc, wvec, avals, ωMvals, γ, Γ, nγ, printout=printout) 

    printout ? output("\n## Diagonalization") : nothing
    printout ? output("BLAS Config: {}", BLAS.get_config()) : nothing
    printout ? output("BLAS Number of Threads: {}", BLAS.get_num_threads()) : nothing
    printout ? output("\nDiagonalization started... ", ending="") : nothing
    t = @elapsed begin
        e, Uix = eigen(H)
    end
    printout ? output("Done in {:10.5f} seconds.", t) : nothing

    ### RETURN SYSTEM OBJECT

    # We assume the photon basis come first, thus phot_range = 1:Np (Np being the number of photon states)
    # which in turn is Np = 2Nc + 1 (Nc is the number of positive momentum photonic states)
    # Similarly, the molecular range starts at Np+1 and ends at Np + Nm, where Nm is the number of molecules
    # That is 2Nc + 2 -> 2Nc + Nm + 1
    # This is the range of pure photonic states in the basis used to construct the Hamiltonian
    rphot = 1:(Ncm)
    # This is the range of pure molecular states in the basis used to construct the Hamiltonian
    rmol = (Ncm+1):(Ncm+Nm)

    # Construct and return the `QuantumWire` system object.
    return LeakyQuantumWire(Uix, e, rphot, rmol, ωMvals, avals, ωc, wvec, γ, Γ)
end
