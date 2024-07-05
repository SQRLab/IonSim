module MolmerSorensen

using QuantumOptics
using IonSim

function construct_MS_chamber(
    # Ion-chain parameters
    trap_frequency,

    # Laser parameters
    intensity,
    wavelength,
    detuning,
    ac_correction,

    # Trap parameters
    B_strength,

    )

    chain = LinearChain(
        ions = [Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")]), Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])], 
        comfrequencies = (x = 3e6, y = 3e6, z = trap_frequency), 
        selectedmodes = (;z = [1],)
    )

    laser1 = Laser(pointing=[(1, 1.), (2, 1.)]) # <-- tuple = (integer: index of ion, float: amplitude at position of that ion (as fraction of predefined intensity))
    laser1.ϵ = x̂
    laser1.k = ẑ
    intensity!(laser1, intensity)
    wavelength!(laser1, wavelength)

    laser2 = Laser(pointing=[(1, 1.), (2, 1.)]) # So here, laser hitting each ion with full intensity
    laser2.ϵ = x̂
    laser2.k = ẑ
    intensity!(laser2, intensity)
    wavelength!(laser2, wavelength)

    d = ac_correction
    Δ_blue = trap_frequency + detuning - d 
    Δ_red = -(trap_frequency + detuning - d )
    laser1.Δ = Δ_blue # blue sideband
    laser2.Δ = Δ_red; # red sideband

    chamber = Chamber(iontrap=chain, B=B_strength, Bhat=(x̂ + ẑ)/√2, lasers=[laser1, laser2]);
    return chamber
end


end