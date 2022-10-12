using PolaritonicSystems
using Test

@testset "Observables" begin

    sys = QuantumWire(
        ΩR = 0.3eV,
        ϵ  = 3.0,
        Nc = 100,
        Nm = 300,
        a  = 10nm,
        σa = 0,
        ωM = 1.0eV,
        σM = 0,
        Ly = 1e7nm,
        Lz = 1e7nm,
        nz = 1,
        ny = 1
    )

    @testset "Wavepacket Properties" begin
        # Create wavepacket and check if it has the right uncertainty
        wvp = create_exciton_wavepacket(1500nm, 80nm, sys)
        @test 79.9 ≤ sqrt(mean_square_disp(wvp, 1500nm, sys)) ≤ 80.1

        # Probability distribution of finding the energy on a given molecule
        pm = get_exciton_prob(wvp, sys)

        # Must recover wavepacket properties
        # Average position
        @test sum(pm .* sys.mol_positions) ≈ 1500.0
        # Variance
        @test sum(pm .* (sys.mol_positions .^ 2)) - sum(pm .* sys.mol_positions)^2 ≈ 6400.0

        # This exciton wavepacket has probability of zero for finding a photon
        @test prob_any_phot(wvp, sys) < 1e-10
        @test prob_any_mol(wvp, sys) ≈ 1.0

        # Initial exciton survival must be one
        @test exciton_survival_prob(sys, 0.0) == 1.0
    end

    # Compare implemented `exciton_survival_prob` to an alternative implementation
    # This implementation uses
    # Pm = (1/Nm) ∑ⱼ | ∑ₓ exp(-iExt/ħ) ⟨j|x⟩⟨x|j⟩ |²
    ħ = 0.0006582119569509067
    @testset "Exciton Survival" begin
        for t = 0.0:0.2:2.0
            out = 0.0
            for j = sys.mol_range
                _c2 = 0.0
                for χ in eachindex(sys.evals)
                    _c2 += exp(-im * sys.evals[χ]*t / ħ) * abs2.(sys.Uix[j, χ])
                end
                out += abs2(_c2)
            end
            @test out / length(sys.mol_positions) ≈ exciton_survival_prob(sys, t)
        end
    end
end
