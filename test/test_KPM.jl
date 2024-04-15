using PolaritonicSystems
using PolaritonicSystems.KPM
using Test

@testset "KPM" begin
    @testset "Type Stability" begin
        for arraytype in ["Dense", "Sparse", "BlockArrowhead"]
            @testset "$arraytype Array" begin
                for T in [Float64, Float32]

                    H = PolaritonicSystems.build_qw_hamiltonian(ΩR=0.1, Nm=100, a=10, Em=2.0, Ecmax=3.5, ny=1, nz=1, Ly=200, Lz=400, ϵ=3, htype=arraytype, dtype=T)
                    @test eltype(H) == complex(T)

                    tH, a, b = KPM.cheb_scale(H)
                    @test eltype(tH) == complex(T)

                    ρ = KPM.DOS(tH, 100, 20)
                    @test eltype(ρ) == T

                    ρloc = KPM.LDOS(tH, 50, 100)
                    @test eltype(ρloc) == T

                    ρrange = KPM.LDOS(tH, 1:20, 100, 20)
                    @test eltype(ρrange) == T
                end
            end
        end
    end
end