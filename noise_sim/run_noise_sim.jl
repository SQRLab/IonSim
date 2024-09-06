# Run this cell to test if the LinearChain object is gonna give us a headache.
using IonSim

@time chain = LinearChain(
        ions=[Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")]), Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])], 
        comfrequencies=(x=3e6,y=3e6,z=2.5e5), selectedmodes=(;z=[1],)
    )
chain = Nothing;

using QuantumOptics
import PyPlot
const plt = PyPlot;
using Random, Distributions
using ProgressBars
using DelimitedFiles
Random.seed!(0)

# Locals
include("./noise_sim.jl")
using .NoiseSim: construct_MS_chamber, get_MS_noise, get_single_qubit_gate_noise

########## Calcium-40 ion ##########

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

chamber_temp = construct_MS_chamber(I=μ_I_MS, ν=μ_ν_MS, ν_target=μ_ν_MS, f_cl=μ_f_cl_MS, ϕ=μ_ϕ_MS, π2_time=π2_TIME_MS, B=B_STRENGTH_MS, ac_correction=AC_CORRECTION)
ket_0_mot = IonSim.modes(chamber_temp)[1][0]
chamber_temp = Nothing

########## Other useful globals ##########
C0 = 2.99792458e8 # Speed of light, m/s
δλ_MAX = 1e-15 # Larger than this and the simulation becomes unstable
N_SHOTS = Int(1e3);


function plot_populations(tout, sol)

    # compute expectation values
    ρ_00 = dm(ket_00 ⊗ ket_0_mot)
    ρ_01 = dm(ket_01 ⊗ ket_0_mot)
    ρ_10 = dm(ket_10 ⊗ ket_0_mot)
    ρ_11 = dm(ket_11 ⊗ ket_0_mot)
    ρ_00_p_i11 = dm(ket_00_p_i11 ⊗ ket_0_mot)
    ρ_00_m_i11 = dm(ket_00_m_i11 ⊗ ket_0_mot)
    ρ_01_p_i10 = dm(ket_01_p_i10 ⊗ ket_0_mot)
    ρ_01_m_i10 = dm(ket_01_m_i10 ⊗ ket_0_mot)

    prob_00 = expect(ρ_00, sol)  # 𝔼(|S⟩|S⟩)
    prob_11 = expect(ρ_11, sol)  # 𝔼(|D⟩|D⟩)
    prob_01 = expect(ρ_01, sol)  # 𝔼(|S⟩|D⟩)
    prob_10 = expect(ρ_10, sol)  # 𝔼(|D⟩|S⟩)
    prob_00_p_i11 = expect(ρ_00_p_i11, sol)  # 𝔼(|S⟩|S⟩ + i|D⟩|D⟩)
    prob_00_m_i11 = expect(ρ_00_m_i11, sol)  # 𝔼(|S⟩|S⟩ - i|D⟩|D⟩)
    prob_01_p_i10 = expect(ρ_01_p_i10, sol)  # 𝔼(|S⟩|D⟩ + i|D⟩|S⟩)
    prob_01_m_i10 = expect(ρ_01_m_i10, sol)  # 𝔼(|S⟩|D⟩ - i|D⟩|S⟩)

    # plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.plot(tout, prob_00, label="00")
    ax1.plot(tout, prob_11, label="11")
    ax1.plot(tout, prob_01, label="01")
    ax1.plot(tout, prob_10, label="10")
    ax1.set_xlim(tout[1], tout[end])
    ax1.set_ylim(0, 1)
    ax1.legend(loc=1)
    ax1.set_xlabel("Time (μs)")
    ax1.set_ylabel("Population")
    ax1.set_title("Computational basis")

    ax2.plot(tout, prob_00_p_i11, label="00 + i11")
    ax2.plot(tout, prob_00_m_i11, label="00 - i11")
    ax2.plot(tout, prob_01_p_i10, label="01 + i10")
    ax2.plot(tout, prob_01_m_i10, label="01 - i10")
    ax2.set_xlim(tout[1], tout[end])
    ax2.set_ylim(0, 1)
    ax2.legend(loc=1)
    ax2.set_xlabel("Time (μs)")
    ax2.set_ylabel("Population")
    ax2.set_title("Bell basis")

    return fig
end

#################################################### Good ####################################################

σ_I_good = 2.5e2 # W/m^2
σ_f_cl_good = 2.5e2 # Hz
σ_ν_good = 2.5e2 # Hz
σ_ϕ_good = 0.005π # rad

### MS gate
println("\n\n################# MS Good ########################")
Normal_I_good_MS = Normal(μ_I_MS, σ_I_good)
Normal_f_cl_good_MS = Normal(μ_f_cl_MS, σ_f_cl_good)
Normal_ν_good_MS = Normal(μ_ν_MS, σ_ν_good)
Normal_ϕ_good_MS = Normal(μ_ϕ_MS, σ_ϕ_good)

fidelities_good_MS, fidelities_el_good_MS, entropies_good_MS = get_MS_noise(I_dist = Normal_I_good_MS, 
                                            f_cl_dist = Normal_f_cl_good_MS, 
                                            ν_dist = Normal_ν_good_MS, 
                                            ϕ_dist = Normal_ϕ_good_MS, 
                                            π2_time = π2_TIME_MS, 
                                            B = B_STRENGTH_MS, 
                                            ac_correction = AC_CORRECTION, 
                                            ν_target = μ_ν_MS, 
                                            n_shots=N_SHOTS)

writedlm("noise_sim_out/fidelities_good_MS_n$N_SHOTS.csv",  fidelities_good_MS, ',')
writedlm("fidelities_el_good_MS_n$N_SHOTS.csv",  fidelities_el_good_MS, ',')
writedlm("entropies_good_MS_n$N_SHOTS.csv",  entropies_good_MS, ',')
println("Min fidelity = $(round(minimum(fidelities_good_MS); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_good_MS); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_good_MS); digits=8))")

### RX
println("\n\n################# RX Good ########################")
Normal_I_good_SINGLE = Normal(μ_I_SINGLE, σ_I_good)
Normal_f_cl_good_SINGLE = Normal(μ_f_cl_SINGLE, σ_f_cl_good)
Normal_ν_good_SINGLE = Normal(μ_ν_SINGLE, σ_ν_good)
Normal_ϕ_good_RX = Normal(0, σ_ϕ_good)

fidelities_good_RX, fidelities_el_good_RX, entropies_good_RX = get_single_qubit_gate_noise(gate_type = "RX",
                                    I_dist = Normal_I_good_SINGLE, # Normal_I_good_SINGLE, 
                                    f_cl_dist = Normal_f_cl_good_SINGLE , # Normal_f_cl_good_SINGLE, 
                                    ν_dist = Normal_ν_good_SINGLE, # Normal_ν_good_SINGLE, 
                                    ϕ_dist = Normal_ϕ_good_RX, #Normal_ϕ_good_RX, 
                                    π2_time = π2_TIME_SINGLE, 
                                    B = B_STRENGTH_SINGLE, 
                                    n_shots=N_SHOTS);

writedlm("fidelities_good_RX_n$N_SHOTS.csv",  fidelities_good_RX, ',')
writedlm("fidelities_el_good_RX_n$N_SHOTS.csv",  fidelities_el_good_RX, ',')
writedlm("entropies_good_RX_n$N_SHOTS.csv",  entropies_good_RX, ',')
println("Min fidelity = $(round(minimum(fidelities_good_RX); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_good_RX); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_good_RX); digits=8))")

### RY
println("\n\n################# RY Good ########################")
Normal_ϕ_good_RY = Normal(π/2, σ_ϕ_good)

fidelities_good_RY, fidelities_el_good_RY, entropies_good_RY = get_single_qubit_gate_noise(gate_type = "RY",
                                    I_dist = Normal_I_good_SINGLE, # Normal_I_good_SINGLE, 
                                    f_cl_dist = Normal_f_cl_good_SINGLE , # Normal_f_cl_good_SINGLE, 
                                    ν_dist = Normal_ν_good_SINGLE, # Normal_ν_good_SINGLE, 
                                    ϕ_dist = Normal_ϕ_good_RY, #Normal_ϕ_good_RX, 
                                    π2_time = π2_TIME_SINGLE, 
                                    B = B_STRENGTH_SINGLE, 
                                    n_shots=N_SHOTS);

writedlm("fidelities_good_RY_n$N_SHOTS.csv",  fidelities_good_RY, ',')
writedlm("fidelities_el_good_RY_n$N_SHOTS.csv",  fidelities_el_good_RY, ',')
writedlm("entropies_good_RY_n$N_SHOTS.csv",  entropies_good_RY, ',')
println("Min fidelity = $(round(minimum(fidelities_good_RY); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_good_RY); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_good_RY); digits=8))")

#################################################### Bad ####################################################

σ_I_bad = 5e2 # W/m^2
σ_f_cl_bad = 5e2 # Hz
σ_ν_bad = 1e3 # Hz
σ_ϕ_bad = 0.01π # rad

### MS gate
println("\n\n################# MS Bad ########################")
Normal_I_bad_MS = Normal(μ_I_MS, σ_I_bad)
Normal_f_cl_bad_MS = Normal(μ_f_cl_MS, σ_f_cl_bad)
Normal_ν_bad_MS = Normal(μ_ν_MS, σ_ν_bad)
Normal_ϕ_bad_MS = Normal(μ_ϕ_MS, σ_ϕ_bad)

fidelities_bad_MS, fidelities_el_bad_MS, entropies_bad_MS = get_MS_noise(I_dist = Normal_I_bad_MS, 
                                            f_cl_dist = Normal_f_cl_bad_MS, 
                                            ν_dist = Normal_ν_bad_MS, 
                                            ϕ_dist = Normal_ϕ_bad_MS, 
                                            π2_time = π2_TIME_MS, 
                                            B = B_STRENGTH_MS, 
                                            ac_correction = AC_CORRECTION, 
                                            ν_target = μ_ν_MS, 
                                            n_shots=N_SHOTS)

writedlm("fidelities_bad_MS_n$N_SHOTS.csv",  fidelities_bad_MS, ',')
writedlm("fidelities_el_bad_MS_n$N_SHOTS.csv",  fidelities_el_bad_MS, ',')
writedlm("entropies_bad_MS_n$N_SHOTS.csv",  entropies_bad_MS, ',')
println("Min fidelity = $(round(minimum(fidelities_bad_MS); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_bad_MS); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_bad_MS); digits=8))")

### RX
println("\n\n################# RX Bad ########################")
Normal_I_bad_SINGLE = Normal(μ_I_SINGLE, σ_I_bad)
Normal_f_cl_bad_SINGLE = Normal(μ_f_cl_SINGLE, σ_f_cl_bad)
Normal_ν_bad_SINGLE = Normal(μ_ν_SINGLE, σ_ν_bad)
Normal_ϕ_bad_RX = Normal(0, σ_ϕ_bad)

fidelities_bad_RX, fidelities_el_bad_RX, entropies_bad_RX = get_single_qubit_gate_noise(gate_type = "RX",
                                    I_dist = Normal_I_bad_SINGLE, # Normal_I_bad_SINGLE, 
                                    f_cl_dist = Normal_f_cl_bad_SINGLE , # Normal_f_cl_bad_SINGLE, 
                                    ν_dist = Normal_ν_bad_SINGLE, # Normal_ν_bad_SINGLE, 
                                    ϕ_dist = Normal_ϕ_bad_RX, #Normal_ϕ_bad_RX, 
                                    π2_time = π2_TIME_SINGLE, 
                                    B = B_STRENGTH_SINGLE, 
                                    n_shots=N_SHOTS);

writedlm("fidelities_bad_RX_n$N_SHOTS.csv",  fidelities_bad_RX, ',')
writedlm("fidelities_el_bad_RX_n$N_SHOTS.csv",  fidelities_el_bad_RX, ',')
writedlm("entropies_bad_RX_n$N_SHOTS.csv",  entropies_bad_RX, ',')
println("Min fidelity = $(round(minimum(fidelities_bad_RX); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_bad_RX); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_bad_RX); digits=8))")

### RY
println("\n\n################# RY Bad ########################")
Normal_ϕ_bad_RY = Normal(π/2, σ_ϕ_bad)

fidelities_bad_RY, fidelities_el_bad_RY, entropies_bad_RY = get_single_qubit_gate_noise(gate_type = "RY",
                                    I_dist = Normal_I_bad_SINGLE, # Normal_I_bad_SINGLE, 
                                    f_cl_dist = Normal_f_cl_bad_SINGLE , # Normal_f_cl_bad_SINGLE, 
                                    ν_dist = Normal_ν_bad_SINGLE, # Normal_ν_bad_SINGLE, 
                                    ϕ_dist = Normal_ϕ_bad_RY, #Normal_ϕ_bad_RX, 
                                    π2_time = π2_TIME_SINGLE, 
                                    B = B_STRENGTH_SINGLE, 
                                    n_shots=N_SHOTS);

writedlm("fidelities_bad_RY_n$N_SHOTS.csv",  fidelities_bad_RY, ',')
writedlm("fidelities_el_bad_RY_n$N_SHOTS.csv",  fidelities_el_bad_RY, ',')
writedlm("entropies_bad_RY_n$N_SHOTS.csv",  entropies_bad_RY, ',')
println("Min fidelity = $(round(minimum(fidelities_bad_RY); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_bad_RY); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_bad_RY); digits=8))")

#################################################### Mid ####################################################
σ_I_mid = 2.5e2 # W/m^2
σ_f_cl_mid = 2.5e2 # Hz
σ_ν_mid = 1e3 # Hz
σ_ϕ_mid = 0.005π # rad

### MS gate
println("\n\n################# MS Mid ########################")
Normal_I_mid_MS = Normal(μ_I_MS, σ_I_mid)
Normal_f_cl_mid_MS = Normal(μ_f_cl_MS, σ_f_cl_mid)
Normal_ν_mid_MS = Normal(μ_ν_MS, σ_ν_mid)
Normal_ϕ_mid_MS = Normal(μ_ϕ_MS, σ_ϕ_mid)

fidelities_mid_MS, fidelities_el_mid_MS, entropies_mid_MS = get_MS_noise(I_dist = Normal_I_mid_MS, 
                                            f_cl_dist = Normal_f_cl_mid_MS, 
                                            ν_dist = Normal_ν_mid_MS, 
                                            ϕ_dist = Normal_ϕ_mid_MS, 
                                            π2_time = π2_TIME_MS, 
                                            B = B_STRENGTH_MS, 
                                            ac_correction = AC_CORRECTION, 
                                            ν_target = μ_ν_MS, 
                                            n_shots=N_SHOTS)

writedlm("fidelities_mid_MS_n$N_SHOTS.csv",  fidelities_mid_MS, ',')
writedlm("fidelities_el_mid_MS_n$N_SHOTS.csv",  fidelities_el_mid_MS, ',')
writedlm("entropies_mid_MS_n$N_SHOTS.csv",  entropies_mid_MS, ',')
println("Min fidelity = $(round(minimum(fidelities_mid_MS); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_mid_MS); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_mid_MS); digits=8))")

### RX
println("\n\n################# RX Mid ########################")
Normal_I_mid_SINGLE = Normal(μ_I_SINGLE, σ_I_mid)
Normal_f_cl_mid_SINGLE = Normal(μ_f_cl_SINGLE, σ_f_cl_mid)
Normal_ν_mid_SINGLE = Normal(μ_ν_SINGLE, σ_ν_mid)
Normal_ϕ_mid_RX = Normal(0, σ_ϕ_mid)

fidelities_mid_RX, fidelities_el_mid_RX, entropies_mid_RX = get_single_qubit_gate_noise(gate_type = "RX",
                                    I_dist = Normal_I_mid_SINGLE, # Normal_I_mid_SINGLE, 
                                    f_cl_dist = Normal_f_cl_mid_SINGLE , # Normal_f_cl_mid_SINGLE, 
                                    ν_dist = Normal_ν_mid_SINGLE, # Normal_ν_mid_SINGLE, 
                                    ϕ_dist = Normal_ϕ_mid_RX, #Normal_ϕ_mid_RX, 
                                    π2_time = π2_TIME_SINGLE, 
                                    B = B_STRENGTH_SINGLE, 
                                    n_shots=N_SHOTS);

writedlm("fidelities_mid_RX_n$N_SHOTS.csv",  fidelities_mid_RX, ',')
writedlm("fidelities_el_mid_RX_n$N_SHOTS.csv",  fidelities_el_mid_RX, ',')
writedlm("entropies_mid_RX_n$N_SHOTS.csv",  entropies_mid_RX, ',')
println("Min fidelity = $(round(minimum(fidelities_mid_RX); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_mid_RX); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_mid_RX); digits=8))")

### RY
println("\n\n################# RY Mid ########################")
Normal_ϕ_mid_RY = Normal(π/2, σ_ϕ_mid)

fidelities_mid_RY, fidelities_el_mid_RY, entropies_mid_RY = get_single_qubit_gate_noise(gate_type = "RY",
                                    I_dist = Normal_I_mid_SINGLE, # Normal_I_mid_SINGLE, 
                                    f_cl_dist = Normal_f_cl_mid_SINGLE , # Normal_f_cl_mid_SINGLE, 
                                    ν_dist = Normal_ν_mid_SINGLE, # Normal_ν_mid_SINGLE, 
                                    ϕ_dist = Normal_ϕ_mid_RY, #Normal_ϕ_mid_RX, 
                                    π2_time = π2_TIME_SINGLE, 
                                    B = B_STRENGTH_SINGLE, 
                                    n_shots=N_SHOTS);

writedlm("fidelities_mid_RY_n$N_SHOTS.csv",  fidelities_mid_RY, ',')
writedlm("fidelities_el_mid_RY_n$N_SHOTS.csv",  fidelities_el_mid_RY, ',')
writedlm("entropies_mid_RY_n$N_SHOTS.csv",  entropies_mid_RY, ',')
println("Min fidelity = $(round(minimum(fidelities_mid_RY); digits=8))")
println("Mean fidelity = $(round(mean(fidelities_mid_RY); digits=8))")
println("Max fidelity = $(round(maximum(fidelities_mid_RY); digits=8))")

println("\n\nAll done!")

