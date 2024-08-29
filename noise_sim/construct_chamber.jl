module MolmerSorensen

using QuantumOptics
using IonSim
function construct_two_ion_chamber(
    I, # Intensity, W/m^2
    ν, # Trap-frequency, Hz
    ν_target, # Target trap-frequency for computing sideband detuning, Hz
    f_cl, # Center-line frequency, Hz
    ϕ, # Relative phase between red and blue sidebands, radian
    ms_π2_time # Time for an MS(π/2) gate, determines the detuning, s
    ;
    B = 6e-4, # Strength of magnetic field, T
    ac_correction = 0 # AC Stark shift correction, Hz
)
    """
    YOU MUST UPDATE ALL THE LASER PARAMETERS IN ORDER TO PERFORM SINGLE-QUBIT GATES.

    By default, this function builds a chamber that is ready to perform a Molmer-Sorensen gate, i.e. the 
    sideband detuning is set such that it couples the internal electronic-states to the motional-modes of 
    the ion-chain, at a strength such that an MS(π/2) gate is performed in the time specified by `ms_π2_time`.

    The single-qubit gates require zero sideband detuning.
    """
    ######### Ion-chain #########
    chain = LinearChain(
        ions = [Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")]), Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])], 
        comfrequencies = (x = 3e6, y = 3e6, z = ν), 
        selectedmodes = (;z = [1],) # <-- Currently only support for axial modes
    )

    ######### Molmer-Sorensen lasers #########

    # Compute the laser wavelength
    C0 = 2.99792458e8
    λ_cl = C0/f_cl

    # Compute the sideband detuning 
    δ = ν_target + 1/ms_π2_time - ac_correction
    Δ_blue = δ
    Δ_red = -δ

    # Compute pointing vectors
    pointing = [(1, 1.), (2, 1.)]
    # This is NOT the same thing as the Poynting vector
    # pointing::tuple = (integer: index of ion, float: amplitude at position of that ion (as a fraction of the intensity))

    laser_red = Laser(λ=λ_cl, I=I, Δ=Δ_red, ϵ=x̂, k=ẑ, ϕ=ϕ, pointing=pointing)
    laser_blue = Laser(λ=λ_cl, I=I, Δ=Δ_blue, ϵ=x̂, k=ẑ, ϕ=0, pointing=pointing)
    
    ######### Chamber #########
    chamber = Chamber(iontrap=chain, B=B, Bhat=(x̂ + ẑ)/√2, lasers=[laser_red, laser_blue]);
    return chamber
end

function MS(chamber, θ, ψ0, ms_π2_time; timescale=1e-6, lamb_dicke_order=1, rwa_cutoff=Inf)

    t_final = ms_π2_time*θ/(π/2)
    t_range = 0:t_final*1e-3:t_final
    h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff);
    tout, sol = timeevolution.schroedinger_dynamic(t_range/timescale, ψ0, h);

    return tout, sol
end

function prep_for_single_qubit_gate(chamber, ion_idx, π_time)
    """
    This function accepts a two-ion chamber that is already primed to perform an MS gate. 

    It updates the parameters of the lasers to perform the RX gate on the ion specified by `ion_idx`. 

    It also returns the original laser parameters so that they can be reset after the single-qubit gate is performed.

    """
    laser_update_idx = ion_idx
    laser_update = chamber.lasers[laser_update_idx]
    laser_ignore = chamber.lasers[laser_update_idx%2 + 1]

    # Collect the original laser parameters
    og_params = Dict(
        "updated_laser_idx" => laser_update_idx,
        "λ" => laser_update.λ,
        "I_updated" => laser_update.I,
        "I_ignored"=> laser_ignore.I,
        "Δ" => laser_update.Δ,
        "ϕ" => laser_update.ϕ,
        "ϵ" => laser_update.ϵ,
        "k" => laser_update.k,
        "pointing" => laser_update.pointing
    ),
        
    # Update the laser parameters
    detuning!(laser_update, 0)

    λ_tr = transitionwavelength(CALCIUM40, ("S", "D"), chamber)
    wavelength!(laser_update, λ_tr)

    intensity!(laser_update, intensity_from_pitime(laser_update, π_time, CALCIUM40, ("S", "D"), chamber))
    intensity!(laser_ignore, 0)

    polarization!(laser_update, x̂)

    wavevector!(laser_update, ẑ)

    pointing!(laser_update, [(ion_idx, 1.), (ion_idx%2+1, 0.)])
    
    return og_params

end

function restore_og_params(chamber, og_params)
    updated_laser_idx = og_params["updated_laser_idx"]
    laser_updated = chamber.lasers[updated_laser_idx]
    laser_ignored = chamber.lasers[updated_laser_idx%2 + 1]

    wavelength!(laser_updated, og_params["λ"])
    intensity!(laser_updated, og_params["I_updated"])
    intensity!(laser_ignored, og_params["I_ignored"])
    detuning!(laser_updated, og_params["Δ"])
    phase!(laser_updated, og_params["ϕ"])
    polarization!(laser_updated, og_params["ϵ"])
    wavevector!(laser_updated, og_params["k"])
    pointing!(laser_updated, og_params["pointing"])

end

function RX(chamber, ion_idx, θ, ψ0, π_time; timescale=1e-6, lamb_dicke_order=1, rwa_cutoff=Inf)

    """
    This function accepts a chamber that is already primed to perform an MS gate. 

    It updates the parameters of the lasers to perform the RX gate on the ion specified by `ion_idx`. 

    Then it performs the specified RX gate on this ion: with the angle `θ`, starting from the initial state `ψ0`.

    Finally, it resets the laser parameters back to their original values.

    """

    # Update the laser parameters
    og_params = prep_for_single_qubit_gate(chamber, ion_idx, π_time) # Update the laser parameters, this happens in-place for the 'chamber' object

    # Set the necessary phase for the RX gate
    laser_update_idx = ion_idx
    laser_update = chamber.lasers[laser_update_idx]
    phase!(laser_update, 0)

    # Perform the RX gate
    t_final = (θ/π)*π_time
    t_range = 0:t_final*1e-3:t_final
    h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff)
    tout, ψt = timeevolution.schroedinger_dynamic(t_range/timescale, ψ0, h)

    # Reset the laser parameters
    restore_og_params(chamber, og_params)

    return tout, ψt

end

function RY(chamber, ion_idx, θ, ψ0, π_time; timescale=1e-6, lamb_dicke_order=1, rwa_cutoff=Inf)

    # Update the laser parameters
    og_params = prep_for_single_qubit_gate(chamber, ion_idx, π_time) # Update the laser parameters, this happens in-place for the 'chamber' object

    # Set the necessary phase for the RY gate
    laser_update_idx = ion_idx
    laser_update = chamber.lasers[laser_update_idx]
    phase!(laser_update, π/2)

    # Perform the RX gate
    t_final = (θ/π)*π_time
    t_range = 0:t_final*1e-3:t_final
    h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff)
    tout, ψt = timeevolution.schroedinger_dynamic(t_range/timescale, ψ0, h)

    # Reset the laser parameters
    restore_og_params(chamber, og_params)

    return tout, ψt

end

function RZ(chamber, ion_idx, θ, ψ0, π_time; timescale=1e-6, lamb_dicke_order=1, rwa_cutoff=Inf)

    # Update the laser parameters
    og_params = prep_for_single_qubit_gate(chamber, ion_idx, π_time) # Update the laser parameters, this happens in-place for the 'chamber' object
    laser_update_idx = ion_idx
    laser_update = chamber.lasers[laser_update_idx]

    ########## RZ(θ) RX(-π/2) RY(θ) RX(π/2) ##########

    # RX(-π/2) = RX(3π/2)
    phase!(laser_update, 0)
    θ1 = 3π/2
    t_final_1 = (θ1/π)*π_time
    t_range_1 = 0:t_final_1*1e-3:t_final_1
    h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff)
    tout_1, ψ_1 = timeevolution.schroedinger_dynamic(t_range_1, ψ0, h)

    # RY(θ)
    phase!(laser_update, π/2)
    t_final_2 = (θ/π)*π_time
    t_range_2 = 0:t_final_2*1e-3:t_final_2
    h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff)
    tout_2, ψ_2 = timeevolution.schroedinger_dynamic(t_range_2, ψ_1[end], h)

    # RX(π/2)
    phase!(laser_update, 0)
    θ3 = π/2
    t_final_3 = (θ3/π)*π_time
    t_range_3 = 0:t_final_3*1e-3:t_final_3
    h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff)
    tout_3, ψ_3 = timeevolution.schroedinger_dynamic(t_range_3, ψ_2[end], h)

    # Reset the laser parameters
    restore_og_params(chamber, og_params)

    tout = vcat(tout_1, tout_2, tout_3)
    ψt = vcat(ψt_1, ψt_2, ψt_3)
    return tout, ψt

end

end