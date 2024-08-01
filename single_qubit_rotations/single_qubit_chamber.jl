module SingleQubitChamber

using QuantumOptics
using IonSim

CALCIUM40 = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
function construct_single_qubit_chamber(
        # Ion-chain parameters
        trap_frequency,

        # Laser parameters
        pi_time;

        # Trap parameters

        # Default laser parameters
        intensity = Nothing,
        wavelength = Nothing,

        # Default trap parameters
        B_strength = Nothing
    )
    
    chain = LinearChain(
        ions = [CALCIUM40], 
        comfrequencies = (x = 3e6, y = 3e6, z = trap_frequency), 
        selectedmodes = (;z = [1],)
    )

    laser = Laser(Δ=0, ϵ = (x̂ - ẑ)/√2, k = (x̂ + ẑ)/√2, ϕ = 0)

    if B_strength == Nothing
        B_strength = 4e-4
    end
    chamber = Chamber(iontrap=chain, B=B_strength, Bhat=ẑ, δB=0, lasers=[laser]);
    
    if wavelength == Nothing
        λ = transitionwavelength(CALCIUM40, ("S", "D"), chamber)
    else
        λ = wavelength
    end
    wavelength!(laser, λ)

    if intensity == Nothing
        I = intensity_from_pitime(laser, pi_time, CALCIUM40, ("S", "D"), chamber)
    else
        I = intensity
    end
    I = intensity_from_pitime(laser, pi_time, CALCIUM40, ("S", "D"), chamber);
    intensity!(laser, I)

    
    return chamber

end

function RX(θ, chamber, ψ0, pi_time, timescale)

    phase!(chamber.lasers[1], 0)

    h = hamiltonian(chamber, timescale=timescale)
    t_final = (θ/π)*pi_time
    t_range = 0:t_final/100:t_final
    tout, ψt = timeevolution.schroedinger_dynamic(t_range/timescale, ψ0, h)

    return tout, ψt
    
end

function RY(θ, chamber, ψ0, pi_time, timescale)

    phase!(chamber.lasers[1], π/2)

    h = hamiltonian(chamber, timescale=timescale)
    t_final = (θ/π)*pi_time
    t_range = 0:t_final/100:t_final
    tout, ψt = timeevolution.schroedinger_dynamic(t_range/timescale, ψ0, h)

    return tout, ψt
    
end

# RZ(θ) = RX(-π/2) RY(θ) RX(π/2)
function RZ(θ, chamber, ψ0, pi_time, timescale)

    tout_1, ψt_1 = RX(3π/2, chamber, ψ0, pi_time, timescale)

    tout_2, ψt_2 = RY(θ, chamber, ψt_1[end], pi_time, timescale)
    tout_2 .+= tout_1[end]

    tout_3, ψt_3 = RX(π/2, chamber, ψt_2[end], pi_time, timescale)
    tout_3 .+= tout_2[end]

    tout = vcat(tout_1, tout_2, tout_3)
    ψt = vcat(ψt_1, ψt_2, ψt_3)

    return tout, ψt
    
end

end