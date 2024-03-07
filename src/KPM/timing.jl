using HDF5

function generate_data(;order=100, Nvals=[100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000])

    path = joinpath(@__DIR__, "O$(order).h5")
    h5write(path, "NVALS", collect(Nvals))

    tfull = zeros(length(Nvals))
    tlocal = zeros(length(Nvals)) 

    NR = 20

    for i in eachindex(Nvals)
        println("N = $(Nvals[i])")
        H = PolaritonicSystems.KPM.build_H(Nvals[i], 1, 3)
        for r in 1:NR
            print(r, " ")
            tfull[i] += @elapsed PolaritonicSystems.KPM.get_DOS(H, order, 20)
            tlocal[i] += @elapsed PolaritonicSystems.KPM.get_LDOS(H, Nvals[i] ÷ 2, order)
            println()
        end
    end

    tfull = tfull ./ NR
    tlocal = tlocal ./ NR

    h5write(path, "DOS", tfull)
    h5write(path, "LDOS", tlocal)
end

function plot_timings()

    dos = h5read(path, "DOS")
    ldos = h5read(path, "LDOS")
    Nvals = h5read(path, "NVALS")

    lf_dos = get_linear_fit(Nvals, dos)
    lf_ldos = get_linear_fit(Nvals, ldos)

    Nforfit = Nvals[1]:50:Nvals[end]

    #scatter(Nvals, dos, label="DOS")
    scatter(Nvals, ldos, label="LDOS")

    #plot!(Nforfit, [lf_dos[1]*n + lf_dos[2] for n in Nforfit], label="LinFit DOS")
    #plot!(Nforfit, [lf_ldos[1]*n + lf_ldos[2] for n in Nforfit], label="Linfit LDOS")

    plot!()
end

function get_linear_fit(x, y::Vector{T}) where T
    X = zeros(typeof(x[1]), length(x), 2)
    X[:,1] .= 1
    X[:,2] .= x
    a,b = X \ y

    return a,b
end
