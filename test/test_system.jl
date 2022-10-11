# Tests for QuantumWire and Hamiltonian construction
using PolaritonicSystems
using Test
using Unitful
using PhysicalConstants

h = PhysicalConstants.CODATA2018.PlanckConstant
c = PhysicalConstants.CODATA2018.SpeedOfLightInVacuum

# Test unit conversion
@testset "Unit conversion" begin
    sys1 = QuantumWire(
        ΩR = 0.3eV,
        ϵ  = 3.0,
        Nc = 100,
        Nm = 100,
        a  = 20nm,
        σa = 0,
        ωM = 1.0eV,
        σM = 0,
        Ly = 1e7nm,
        Lz = 1e7nm,
        nz = 1,
        ny = 1
    )

    sys2 = QuantumWire(
        ΩR = 241.9663287u"mm^-1",
        ϵ  = 3.0,
        Nc = 100,
        Nm = 100,
        a  = 200u"Å",
        σa = 0,
        ωM = 8065.54429u"cm^-1",
        σM = 0,
        Ly = 1e4u"μm",
        Lz = 10u"mm",
        nz = 1,
        ny = 1
    )

    sys3 = QuantumWire(
        ΩR = 455.7802542u"ps^-1",
        ϵ  = 3.0,
        Nc = 100,
        Nm = 100,
        a  = 20000u"pm",
        σa = 0,
        ωM = 1.5192675143e6u"ns^-1",
        σM = 0,
        Ly = 1e4u"μm",
        Lz = 0.01u"m",
        nz = 1,
        ny = 1
    )

    @test sys1.evals ≈ sys2.evals
    @test sys1.evals ≈ sys3.evals
    @test sys2.evals ≈ sys3.evals
end


@testset "Tavis-Cuming Limit" begin

    # Fix all molecules at x = 0
    # No disorder
    # Adjust length of the cavity such that we have one resonant mode
    ΩR = 0.3eV
    ωM = 1.0eV
    Lyz = h * c / (√2 * ωM)
    sys = QuantumWire(
        ΩR = ΩR,
        ϵ  = 1.0,
        Nc = 0,
        Nm = 100,
        a  = 0nm,
        σa = 0,
        ωM = ωM,
        σM = 0,
        Ly = Lyz,
        Lz = Lyz,
        nz = 1,
        ny = 1
    )

    # Test if the splitting is equal the Rabi Splitting
    @test sys.evals[end] - sys.evals[1] ≈ ustrip(u"eV", ΩR) 
end
