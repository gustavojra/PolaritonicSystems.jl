using Distributions

struct Dipole{T}
    Em::T
    coord::Tuple{T, T, T}
    μ::Tuple{T,T,T}
end

struct PhotonMode{T}
    Ec::T
    k::Tuple{T, T, T}
end

function get_planarcavity_modes(Ncx, Ncy, Lx, Ly, nz, Lz, ϵ; T = Float64,
    c = ustrip(u"nm/ps", CODATA2018.SpeedOfLightInVacuum), ħ = ustrip(u"eV*ps", CODATA2018.PlanckConstant)/2π)

    # Compute the k component along z (short dimension)
    kz = nz*π/Lz

    NT = (2*Ncx + 1) * (2*Ncy + 1)

    out = Array{PhotonMode{T}}(undef, NT)

    i = 1
    for mx in -Ncx:1:Ncx
        qx = 2π * mx / Lx

        for my in -Ncy:1:Ncy
            qy = 2π * my / Ly

            # Wave vector, including z componenet
            k = (qx, qy, kz)

            # Energy
            E = ħ*c*norm(k)/ϵ

            out[i] = PhotonMode(E, k)
            i += 1
        end
    end

    return out
end

function get_planarcavity_molecules(Nm1, Nm2, z, Em, σm, a, σa; distype="normal", T = Float64)

    # Initialize output - array of dipole objects
    out = Array{Dipole{T}}(undef, Nm1*Nm2)

    dist = nothing
    if lowercase(string(distype)) == "normal"
        dist = Normal(Em, σm)
    elseif lowercase(string(distype)) == "uniform"
        dist = Uniform(Em, σm)
    else
        error("Distribution option $dist invalid for static disorder.")
    end

    x = zero(T)
    for i in 1:Nm1
        x += a
        y = zero(T)
        for j in 1:Nm2
            # Get excitation energy
            Em = rand(dist)

            att = 1
            if Em < 0
                if att > 10
                    error("Cannot generate distribution with only positive energies for Em = $Em and σm = $σm")
                end
                Em = rand(dist)
                att += 1
            end

            # Get coordinates
            y += a
            coord = (x + rand(Normal(0,σa)), y + rand(Normal(0,σa)), T(z))

            # Get dipole direction
            θ = rand(Uniform(0, 2π))
            ϕ = rand(Uniform(0, π))
            μ = (sin(ϕ)*cos(θ), sin(ϕ)*sin(θ), cos(ϕ))

            # Create and store dipole object
            out[j + Nm2*(i-1)] = Dipole(Em, coord, μ)
        end
    end

    return out
end

function build_pc_hamiltonian(dipoles::Vector{Dipole{T}}, modes::Vector{PhotonMode{T}}, ΩR) where T

    Nm = length(dipoles)
    NT = Nm + 2*length(modes)

    TErange = 1:length(modes)
    TMrange = TErange .+ length(modes)
    drange = (2*length(modes)+1):NT

    H = zeros(complex(T), NT, NT)

    # Populate diagonals
    for (i,j) in zip(TErange, TMrange)
        H[i,i] = modes[i].Ec
        H[j,j] = modes[i].Ec
    end

    for (i,d) in zip(drange, dipoles)
        H[i,i] = d.Em
    end

    # Get interaction between dipoles and TE modes
    εz = [0.0, 0.0, 1.0]
    for (i, d) = zip(drange, dipoles)
        for (j, TE) = zip(TErange, modes)

            kz = TE.k[3]
            z = d.coord[3]
            q = [TE.k[1], TE.k[2], 0.0]
            εq = q ./ norm(q)

            # Prefactor, including Rabi splitting
            pf = -im * 0.5 * ΩR * √(d.Em / (Nm*TE.Ec)) 

            # TE Spatial profile - except exp part
            f = sin(kz * z) * cross(εz, εq) 

            # Interaction term assembled
            H[i,j] = pf * dot(d.μ, f)  * exp(im * dot(q, d.coord))
            H[j,i] = adjoint(H[i,j])
        end
    end

    # Get interaction between dipoles and TM modes
    εz = [0.0, 0.0, 1.0]
    for (i, d) = zip(drange, dipoles)
        for (j, TM) = zip(TMrange, modes)

            kz = TM.k[3]
            z = d.coord[3]
            q = [TM.k[1], TM.k[2], 0.0]
            εq = q ./ norm(q)

            # Prefactor, including Rabi splitting
            pf = -im * 0.5 * ΩR * √(d.Em / (Nm*TM.Ec)) 

            # TM Spatial profile - except exp part
            f = (norm(TM.k)/kz)^2 * (sin(kz * z) .* εq - (im * cos(kz * z) / kz) .* εz)

            # Interaction term assembled
            H[i,j] = pf * dot(d.μ, f)  * exp(im * dot(q, d.coord))
            H[j,i] = adjoint(H[i,j])
        end
    end

    return H
end

function build_bah_pc_hamiltonian(dipoles::Vector{Dipole{T}}, modes::Vector{PhotonMode{T}}, ΩR) where T

    Nm = length(dipoles)
    NTE = length(modes)
    NT = Nm + 2*NTE

    d = zeros(complex(T), NT)
    X = zeros(complex(T), Nm, 2*length(modes))

    # Populate diagonals
    for i in eachindex(modes)
        d[i] = modes[i].Ec
        d[i+NTE] = modes[i].Ec
    end

    for i in eachindex(dipoles)
        d[i + 2*NTE] = dipoles[i].Em
    end

    εz = [0.0, 0.0, 1.0]
    for (i, dip) = enumerate(dipoles)

        # Get interaction between dipoles and TE modes
        for (j, M) = enumerate(modes)

            kz = M.k[3]
            z = dip.coord[3]
            q = [M.k[1], M.k[2], 0.0]
            εq = norm(q) > 0.0 ? q ./ norm(q) : [1.0, 0.0, 0.0]

            # Prefactor, including Rabi splitting
            pf = -im * 0.5 * ΩR * √(dip.Em / (Nm*M.Ec)) 

            # TE Spatial profile - except exp part
            fTE = sin(kz * z) * cross(εz, εq) 

            # Interaction term assembled
            X[i,j] = pf * dot(dip.μ, fTE)  * exp(im * dot(q, dip.coord))
            #X[i,j] = pf * norm(fTE)  * exp(im * dot(q, dip.coord))

            # TM Spatial profile - except exp part
            fTM = (1/norm(M.k)) * (kz * sin(kz * z) .* εq - (im * norm(q) * cos(kz * z) / kz) .* εz)

            # Interaction term assembled
            X[i,j+NTE] = pf * dot(dip.μ, fTM)  * exp(im * dot(q, dip.coord))
            #X[i,j+NTE] = pf * norm(fTM)  * exp(im * dot(q, dip.coord))
        end
    end

    return SymBlockArrowHead(d, X, Nm+2*NTE, 1:(2*NTE), (2*NTE+1):(NT))
end






