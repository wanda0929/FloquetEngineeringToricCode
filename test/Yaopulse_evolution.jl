using Yao
using CairoMakie
using LinearAlgebra
using FloquetEngineeringToricCode
using DifferentialEquations
using CairoMakie: PolarAxis
#the first simulation
reg1, reg2, res1, res2, evo_1 = sim1(interaction=false)
fig = plot_results(res1, res2, pulses=[4])

#extract the effective hamiltonian
omega = 1
tau = 2*pi/omega
#check if the evolution operator is correct
#evo_1 = operator_evolve(pulse_hamiltonian_interaction_1(1, 1.3, 2.6, 0.35), 2π, 10000)
#evo_1 = operator_evolve(pulse_hamiltonian_interaction(1, 12.0, 10.0, 1.3, 2.6, 0.35), 2π, 10000)

rr = statevec(zero_state(4))
rr1 = evo_1 * rr
statevec(reg1)
diff = real.(rr1 - statevec(reg1))  # 提取实数部分
println("the difference between the two is: ", diff)
println("Maximum absolute difference: ", maximum(abs.(diff)))
XZZX = XZZX_operator(tau, π/8/tau, 10000)
# 将矩阵转换为 AbstractBlock
evo_block = matblock(evo_1)
XZZX_block = matblock(XZZX)
operator_fidelity(evo_block, XZZX_block)

#extract the effective hamiltonian
omega = 1
tau = 2*pi/omega
evo = log(evo_1) * im / (tau)

# 绘制 Pauli 矩阵迹的极坐标图
fig_polar = plot_pauli_traces(evo, title="Pauli Matrix Traces (First Simulation)")


#the second simulation
reg1, res1, res2, evo_2 = sim2(; interaction=false)
# Plot the results from sim2
fig = plot_results(res1, res2, pulses=[4])

#extract the effective hamiltonian
evo_2 = operator_evolve(pulse_h(1, 14.0, 11.0, 1.2, 2.6, 0.054), 2π, 10000)
omega = 1
tau = 2*pi/omega
evoo = log(evo_2) * im / (tau)

#check if the evolution operator is correct
rr = statevec(zero_state(4))
rr1 = evo_2 * rr
statevec(reg1)
diff = real.(rr1 - statevec(reg1))  # 提取实数部分
println("the difference between the two is: ", diff)
println("Maximum absolute difference: ", maximum(abs.(diff)))

# 绘制第二个模拟的 Pauli 矩阵迹的极坐标图
fig_polar2 = plot_pauli_traces(evoo, title="Pauli Matrix Traces (Second Simulation)")
display(fig_polar2)