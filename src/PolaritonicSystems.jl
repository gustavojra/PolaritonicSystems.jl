module PolaritonicSystems
using PhysicalConstants: CODATA2018
using Formatting
using Unitful

export eV, nm
export output
export QuantumWire, time_propagate, create_exciton_wavepacket, mean_square_disp, mean_disp
export time_propagate!, eigenvector_spread
export prob_any_mol, prob_any_phot, get_exciton_prob, exciton_survival_prob
export average_energy, average_square_energy, energy_uncertainty, energy_filter!
export LeakyQuantumWire, create_photon_wavepacket

const ħ = ustrip(u"eV*ps", CODATA2018.PlanckConstant) / 2π
const imħ = -im/ħ

function _print_inputs(inputs...)

    output("           ---Quantum System Input Parameters---       ")
    output("\n              # Light-Matter Interaction #             ")
    output("Rabi Frequency (ΩR):                                   {}", inputs[1])

    output("\n                # Molecular Properties #               ")
    output("Number of molecules (Nm):                              {}", inputs[2])
    output("Average intermolecular distance (a):                   {}", inputs[3])
    output("Intermolecular distance standard deviation (σa):       {}", inputs[4])
    output("Average molecular excitation energy (ωM):              {}", inputs[5])
    output("Molecular excitation energy standard deviation (σM):   {}", inputs[6])

    output("\n                # Cavity/Wire Properties #               ")
    output("Number of positive cavity modes (Nc):                  {}", inputs[7])
    output("Cavity y dimension size (Ly):                          {}", inputs[8])
    output("Cavity z dimension size (Lz):                          {}", inputs[9])
    output("y-momentum quantum number (ny):                        {}", inputs[10])
    output("z-momentum quantum number (nz):                        {}", inputs[11])
    output("Index of refraction squared (ϵ):                       {}", inputs[12])
end

function output(str, x...; ending="\n")
    f = format(str*ending, x...)

    open("log.out", "a") do io
        write(io, f)
        flush(io)
    end
end

# Hamiltonian Constructors
include("Hamiltonian.jl")

# Quantum System Objects
include("QuantumWire.jl")
include("PlanarCavity.jl")

# States and time evolution
include("States.jl")

# Observables (measurements) and time propagation is implemented here
include("Observables.jl")

include("KPM/KPM.jl")

end # module
