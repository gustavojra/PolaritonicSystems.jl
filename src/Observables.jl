export mean_r2

"""
    mean_square_disp

Computes the average squared displacement (with respect to an initial position x0): ⟨(x - x0)²⟩

# Arguments

| Argument |       Type      |                                       Description                                                                  |
|:--------:|:----------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector          | Vector representing the state to be analyzed.                                                                      |
| x0       | Number/Quantity | Position for which the displacement is computed against. If a number, unit is taken as `nm`.                       |
| system   | QuantumWire     | QuantumWire object - constaints all static information of the system.                                              |
"""
function mean_square_disp(state::Vector, x0::Number, system::QuantumWire)

    # Transform to localized (uncoupled) state
    unc_state = system.Uix * state

    # Loop through molecules
    x2 = 0.0
    n = 1
    Pmol = 0.0
    for i = system.mol_range
        # Contribution of each molecule is P * (x-x0)^2, where P is the probability of the state
        # being localized at that molecule and x is its position
        x2 += abs2(unc_state[i]) * (system.mol_positions[n] - x0)^2

        # Probabilty of finding a molecule excited
        Pmol += abs2(unc_state[i])
        n += 1
    end

    # Renormalize (conditional probability) with the probability of finding the molecular exciton
    return x2 / Pmol
end

# Support for units
function mean_square_disp(state::Vector, x0::Quantity, system::QuantumWire)
    mean_square_disp(state, ustrip(u"nm", x0), system)
end

"""
    mean_disp

Computes the average displacement (with respect to an initial position x0): ⟨(x - x0)⟩

# Arguments

| Argument |       Type      |                                       Description                                                                  |
|:--------:|:----------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector          | Vector representing the state to be analyzed.                                                                      |
| x0       | Number/Quantity | Position for which the displacement is computed against. If a number, unit is taken as `nm`.                       |
| system   | QuantumWire     | QuantumWire object - constaints all static information of the system.                                              |
"""
function mean_disp(state::Vector, x0::Number, system::QuantumWire)
    # Transform to localized (uncoupled) state
    unc_state = system.Uix * state

    # Loop through molecules
    x = 0.0
    n = 1
    Pmol = 0.0
    for i = system.mol_range
        # Contribution of each molecule is P * (x-x0), where P is the probability of the state
        # being localized at that molecule and x is its position
        x += abs2(unc_state[i]) * (system.mol_positions[n] - x0)

        # Probabilty of finding a molecule excited
        Pmol += abs2(unc_state[i])
        n += 1
    end

    # Renormalize (conditional probability) with the probability of finding the molecular exciton
    return x / Pmol
end

# Support for units
function mean_disp(state::Vector, x0::Quantity, system::QuantumWire)
    mean_disp(state, ustrip(u"nm", x0), system)
end

"""
    prob_any_phot

For a given state vector, computes the probability of finding a photon in any mode.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function prob_any_phot(state::Vector, system::QuantumWire)
    locstate = system.Uix * state 
    return abs2.(locstate[system.phot_range]) |> sum
end

"""
    prob_any_mol

For a given state vector, computes the probability of finding an exciton on any molecule.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function prob_any_mol(state::Vector, system::QuantumWire)
    locstate = system.Uix * state 
    return abs2.(locstate[system.mol_range]) |> sum
end

"""
    get_exciton_prob

For a given state vector, returns a renormalized vector with probabilities of finding the excitation energy on each molecule.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function get_exciton_prob(state::Vector, system::QuantumWire)
    loc = system.Uix * state
    loc = abs2.(loc[system.mol_range]) 

    # Renormalize 
    return loc ./ sum(loc)
end

function mean_r2(wvp, dipoles::Vector{Dipole{T}}, r0, drange) where T

    P = sum(abs2.(wvp[drange]))

    out = 0.0
    for (i,d) in zip(drange, dipoles)
        r = sum((r0 .- d.coord) .^ 2)
        out += abs2(wvp[i]) * r
    end

    return out, P
end

"""
    excition_survival_prob

For a given system, compute the mean probabiility that an initially excited molecule at time t = 0 will remain excited at time t. 
Pₘ(t) = ∑ₙ |⟨1ₙ|exp(-iĤt/ħ)|1ₙ⟩|² / Nm

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
| t        | Number/Quantity | Time for which the system must be evolved into. If a number, unit is taken as `ps`.                            |
"""
function exciton_survival_prob(sys::QuantumWire, t::Quantity)
    exciton_survival_prob(sys, ustrip(u"ps", t))
end

function exciton_survival_prob(sys::QuantumWire, t::Number)
    out = 0.0
    expv = [exp(imħ*E*t) for E in sys.evals]

    for i = sys.mol_range
        dot = 0.0
        # Compute |⟨initial|final⟩|²
        for j = eachindex(sys.evals)
            c = sys.Uix[i,j]
            dot += adjoint(c) * expv[j] * c
        end
        out += abs2(dot)
    end

    # Divide it by Nm
    return out / length(sys.mol_range)
end

"""
    average_energy

For a given state vector, returns the average energy  ⟨E⟩ = ⟨ψ|Ĥ|ψ⟩.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function average_energy(state::Vector, sys::QuantumWire)
    return dot(state, sys.evals .* state) |> real
end

"""
    average_square_energy

For a given state vector, returns the average squared energy  ⟨E²⟩ = ⟨ψ|Ĥ²|ψ⟩.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function average_square_energy(state::Vector, sys::QuantumWire)
    return dot(state, (sys.evals.^2) .* state) |> real
end

"""
    average_square_energy

For a given state vector, returns the energy uncertainty σₑ = √[ ⟨E²⟩ - ⟨E⟩² ] 

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function energy_uncertainty(state::Vector, sys::QuantumWire)
    E2 = average_square_energy(state, sys)
    E = average_energy(state, sys)

    # Abs is needed because sometimes E2 is just numerically smaller than E² and we get a domain error for sqrt
    return √(abs(E2 - E^2))
end

function eigenvector_spread(sys::QuantumWire)
    out = zeros(length(sys.evals))

    for χ in eachindex(out)
        r = 0.0
        r2 = 0.0
        for (n,i) in enumerate(sys.mol_range)
            r  += abs2(sys.Uix[i, χ]) * sys.mol_positions[n]
            r2 += abs2(sys.Uix[i, χ]) * sys.mol_positions[n]^2
        end

        out[χ] = sqrt(r2 - r^2)
    end

    return out
end

function DOS(sys::QuantumWire, δE)
    return DOS(sys.evals, δE)
end

function DOS(H, δE)
    e, _ = eigen(H)
    return DOS(e, δE)
end

function DOS(evals::Vector, δE)

    # Find boundaries for DOS evaluation
    Emin = minimum(evals)
    Emax = maximum(evals) + 0.01 # Add a little shift here to ensure the largest eigenvalue is included in the last bin

    # Range containing lower and upper limits for each energy bin
    Elims = range(start=Emin, stop=Emax, step=δE)

    # Initialize dos array 
    ρ = zeros(length(Elims)-1)
    Ebars = zeros(length(Elims)-1)

    # Loop through bin limits
    for i in eachindex(Elims)

        # Skip the last one (since there would be no i+1 entry). 
        if i == length(Elims)
            break
        end

        # Store the verage value of the bin
        Ebars[i] = 0.5*(Elims[i] + Elims[i+1])

        # Count the number of states within the limits
        ρ[i] = count(x-> x ≥ Elims[i] && x < Elims[i+1], evals)
    end

    # Normalization factor such that the summed area of the histogram yields 1. 
    N = 1 / (δE*sum(ρ))

    # Return bin values and heights. 
    return Ebars, N .* ρ
end

LDOS(sys::QuantumWire, x, δE) = LDOS(sys.evals, sys.Uix, x, δE)
LDOS(sys::QuantumWire, δE) = LDOS(sys.evals, sys.Uix, δE)

function LDOS(evals::Vector, U::Matrix, x, δE)

    Emin = minimum(evals)
    Emax = maximum(evals) + 0.01

    Elims = range(start=Emin, stop=Emax, step=δE)

    Ebars = zeros(length(Elims)-1)
    ρ = zeros(length(Ebars), length(evals))

    # Intermediate array
    ldos = similar(Ebars)

    for i in eachindex(Elims)

        # Skip the last one (since there would be no i+1 entry). 
        if i == length(Elims)
            break
        end

        # Store the verage value of the bin
        Ebars[i] = 0.5*(Elims[i] + Elims[i+1])

        # Find the index of all eigenvalues within the limits
        m = findall(x-> x ≥ Elims[i] && x < Elims[i+1], evals)

        if isempty(m)
            continue
        else
            ldos[i] = sum(abs2.(U[x, m]))
        end
    end

    # Normalize each LDOS
    ldos = ldos / (length(evals) * δE * sum(ldos))

    return Ebars, ρ
end

function LDOS(evals::Vector, U::Matrix, δE)

    Emin = minimum(evals)
    Emax = maximum(evals) + 0.01

    Elims = range(start=Emin, stop=Emax, step=δE)

    Ebars = zeros(length(Elims)-1)
    ρ = zeros(length(Ebars), length(evals))

    # Intermediate array
    ldos = similar(Ebars)

    for x in eachindex(evals)

        ldos .= 0.0
        for i in eachindex(Elims)

            # Skip the last one (since there would be no i+1 entry). 
            if i == length(Elims)
                break
            end

            # Store the verage value of the bin
            Ebars[i] = 0.5*(Elims[i] + Elims[i+1])

            # Find the index of all eigenvalues within the limits
            m = findall(x-> x ≥ Elims[i] && x < Elims[i+1], evals)

            if isempty(m)
                continue
            else
                ldos[i] = sum(abs2.(U[x, m]))
            end
        end

        # Normalize each LDOS
        ρ[:,x] = ldos / (length(evals) * δE * sum(ldos))
    end

    return Ebars, ρ
end