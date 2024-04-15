using PolaritonicSystems
using Test

@testset "PolaritonicSystems" begin
    include("test_system.jl")
    include("test_observables.jl")
    include("test_KPM.jl")
end
