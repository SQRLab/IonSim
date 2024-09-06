
module NoiseSim

using QuantumOptics
using IonSim

C0 = 2.99792458e8 # Speed of light, m/s
δλ_MAX = 1e-15 # Larger than this and the simulation becomes unstable
CALCIUM40 = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])

########## Some calibrated values for the Molmer-sorensen gate ##########

μ_I_MS = 94326.65907221894 # laser-intensity, W/cm^2
μ_ν_MS = 2.5e5 # trap-frequency, Hz
μ_f_cl_MS = 4.111550352057269e14 # laser-frequency, Hz
μ_ϕ_MS = 0.0 # relative phase between red and blue sidebands 
π2_TIME_MS = 1e-4 # gate-time for an MS(π/2) gate
AC_CORRECTION = 0.0 # ac-stark shift correction, Hz
B_STRENGTH_MS = 6e-4 # magnetic field strength, T

########## Some calibrated values for the single-qubit gates ##########

μ_I_SINGLE = 4.123759471808808e6
μ_ν_SINGLE = 1e6 # Although it shouldn't really matter
μ_f_cl_SINGLE = 4.111550340833542e14
μ_ϕ_SINGLE = 0.0
π2_TIME_SINGLE = 1e-6
B_STRENGTH_SINGLE = 4e-4

########## Some useful quantum states ##########
ket_0 = CALCIUM40["S"]
ket_1 = CALCIUM40["D"]

# Computational basis
ket_00 = CALCIUM40["S"] ⊗ CALCIUM40["S"]
ket_01 = CALCIUM40["S"] ⊗ CALCIUM40["D"]
ket_10 = CALCIUM40["D"] ⊗ CALCIUM40["S"]
ket_11 = CALCIUM40["D"] ⊗ CALCIUM40["D"]
ρ_00 = dm(ket_00)
ρ_01 = dm(ket_01)
ρ_10 = dm(ket_10)
ρ_11 = dm(ket_11)

# Bell basis
ket_00_m_i11 = (ket_00 - 1im*ket_11)/√2
ket_00_p_i11 = (ket_00 + 1im*ket_11)/√2
ket_01_m_i10 = (ket_01 - 1im*ket_10)/√2
ket_01_p_i10 = (ket_01 + 1im*ket_10)/√2
ρ_00_m_i11 = dm(ket_00_m_i11)
ρ_00_p_i11 = dm(ket_00_p_i11)
ρ_01_m_i10 = dm(ket_01_m_i10)
ρ_01_p_i10 = dm(ket_01_p_i10);


function construct_MS_chamber(;
    I = Nothing, # Intensity, W/m^2
    ν = Nothing, # Trap-frequency, Hz
    ν_target = Nothing, # Target trap-frequency for computing sideband detuning, Hz
    f_cl = Nothing, # Center-line frequency, Hz
    ϕ = Nothing, # Relative phase between red and blue sidebands, radian
    π2_time = Nothing, # Time for an MS(π/2) gate, determines the detuning, s
    B = Nothing, # Strength of magnetic field, T
    ac_correction = Nothing # AC Stark shift correction, Hz
)

    ######### Ion-chain #########
    chain = LinearChain(
        ions = [Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")]), Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])], 
        comfrequencies = (x = 3e6, y = 3e6, z = ν), 
        selectedmodes = (;z = [1],) # <-- Currently only support for axial modes
    )

    ######### Molmer-Sorensen lasers #########

    # Compute the laser wavelength
    λ_cl = C0/f_cl

    # Compute the sideband detuning 
    δ = ν_target + 1/π2_time - ac_correction
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

chamber_temp = construct_MS_chamber(I=μ_I_MS, ν=μ_ν_MS, ν_target=μ_ν_MS, f_cl=μ_f_cl_MS, ϕ=μ_ϕ_MS, π2_time=π2_TIME_MS, B=B_STRENGTH_MS, ac_correction=AC_CORRECTION)
ket_0_mot = IonSim.modes(chamber_temp)[1][0]
chamber_temp = Nothing

function get_MS_noise(;
    I_dist = Nothing, 
    f_cl_dist = Nothing, 
    ν_dist = Nothing, 
    ϕ_dist = Nothing, 
    π2_time = Nothing, 
    B = Nothing, 
    ac_correction = Nothing, 
    ν_target = Nothing,
    n_shots=1000, timescale=1e-6, lamb_dicke_order=1, rwa_cutoff=Inf)
    # Initial state
    ψ0 = ket_00 ⊗ ket_0_mot

    # Target state
    ket_el_target = ket_00_p_i11
    ρ_target = dm(ket_el_target ⊗ ket_0_mot)
    ρ_el_target = ptrace(ρ_target, 3)

    # Simulation time
    t_range = 0:π2_time*1e-3:π2_time

    # Collect data
    fidelities = []
    fidelities_el = []
    entropies = []
    n_collected = 0
    while n_collected < n_shots

        I = rand(I_dist)
        # Intensity can't be negative
        if I<0
            println("I<0")
            continue
        end

        f_cl = rand(f_cl_dist)
        λ_cl = C0/f_cl
        # If we're too far from the target wavelength, simulation becomes unstable
        δλ = abs(λ_cl - C0/μ_f_cl_MS)
        if δλ > 1e-15
            println("δλ > 1e-15")
            continue
        end

        ν = rand(ν_dist)
        ϕ = rand(ϕ_dist)

        # Construct chamber with given parameters
        try
            chamber = construct_MS_chamber(I=I, ν=ν, ν_target=ν_target, f_cl=f_cl, ϕ=ϕ, π2_time=π2_time, B=B, ac_correction=ac_correction)

            # Perform gate
            h = hamiltonian(chamber, timescale=timescale, lamb_dicke_order=lamb_dicke_order, rwa_cutoff=rwa_cutoff);
            tout, sol = timeevolution.schroedinger_dynamic(t_range/timescale, ψ0, h);

            # Evaluate performance
            ρ = dm(sol[end])
            ρ_el = ptrace(ρ, 3)
            fid = real(fidelity(ρ, ρ_target))
            fid_el = real(fidelity(ρ_el, ρ_el_target))

            append!(fidelities, fid)
            append!(fidelities_el, fid_el)
            append!(entropies, entropy_vn(ρ_el))
        catch AssertionError # solver failed to find normal modes compatible with characteristic frequencies
            println("AssertionError")
            continue
        end

        n_collected = length(fidelities)
        if n_collected%100 == 0
            println("n_collected = ", n_collected)
        end

    end

    return fidelities, fidelities_el, entropies

end

function construct_single_qubit_chamber(;
    I = Nothing,
    ν = Nothing, # Trap-frequency, Hz
    f_cl = Nothing,
    ϕ = Nothing, # Phase, rad
    π2_time = Nothing, # Time for a π-pulse, s
    B = Nothing # Strength of magnetic field, T
)

    # Ion-chain
    chain = LinearChain(
        ions = [CALCIUM40], 
        comfrequencies = (x = 3e6, y = 3e6, z = ν), 
        selectedmodes = (;z = [1],)
    )

    # Laser
    λ = C0/f_cl
    laser = Laser(I=I, λ=λ, Δ=0, ϵ = (x̂ - ẑ)/√2, k = (x̂ + ẑ)/√2, ϕ = ϕ)

    # Ion trap
    chamber = Chamber(iontrap=chain, B=B, Bhat=ẑ, δB=0, lasers=[laser]);

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

function get_single_qubit_gate_noise(;
    gate_type = Nothing,
    I_dist = Nothing, 
    f_cl_dist = Nothing, 
    ν_dist = Nothing, 
    ϕ_dist = Nothing, 
    π2_time = Nothing, 
    B = Nothing, 
    n_shots=1000, timescale=1e-6, lamb_dicke_order=1, rwa_cutoff=Inf)

    if gate_type=="RX"
        gate_func = RX
        ket_el_target = (ket_0 - 1im*ket_1)/√2
    elseif gate_type=="RY"
        gate_func = RY
        ket_el_target = (ket_0 - ket_1)/√2
    else
        println("Invalid gate type")
        return
    end

    # Initial state
    ψ0 = ket_0 ⊗ ket_0_mot

    # Target state
    ρ_target_el = dm(ket_el_target)
    ρ_target = dm(ket_el_target ⊗ ket_0_mot)

    # Collect data
    fidelities = []
    fidelities_el = []
    entropies = []
    n_collected = 0
    while n_collected < n_shots
        I = rand(I_dist)
        # Intensity can't be negative
        if I<0
            continue
        end

        f_cl = rand(f_cl_dist)
        λ_cl = C0/f_cl
        # If we're too far from the target wavelength, simulation becomes unstable
        δλ = abs(λ_cl - C0/μ_f_cl_SINGLE)
        if δλ > 1e-15
            continue
        end

        ν = rand(ν_dist)
        ϕ = rand(ϕ_dist)

        try
            # Construct chamber with given parameters
            chamber = construct_single_qubit_chamber(I=I, ν=ν, f_cl=f_cl, ϕ=ϕ, π2_time=π2_time, B=B)

            # Perform gate
            tout, sol = gate_func(π/2, chamber, ψ0, 2*π2_TIME_SINGLE, 1e-6)

            # Evaluate performance
            ρ = dm(sol[end])
            ρ_el = ptrace(ρ, 2)
            fid = real(fidelity(ρ, ρ_target))
            fid_el = real(fidelity(ρ_el, ρ_target_el))
            ee = entropy_vn(ρ_el)

            append!(fidelities, fid)
            append!(fidelities_el, fid_el)
            append!(entropies, ee)
            
        catch AssertionError # solver failed to find normal modes compatible with characteristic frequencies
            continue
        end

        n_collected = length(fidelities)
        
    end

    return fidelities, fidelities_el, entropies

end

end