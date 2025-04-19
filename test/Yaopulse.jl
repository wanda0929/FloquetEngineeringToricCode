using Test
using Yao
using CairoMakie
using LinearAlgebra
using FloquetEngineeringToricCode
@testset "operator_evolve" begin
    evo_1, final_state = operator_evolve(2π/10000, 1.0, 12.0, 10.0, 1.3, 2.6, 0.35, 1)
    reg = zero_state(4)
    ω = 1.0
    τ = 2π/ω
    hamiltonian8 = pulse_hamiltonian(ω, 12.0ω , 10.0ω , 1.3ω , 2.6ω , 0.35ω)
    reg1, res1, reg2, res2, evolution_operator = evolve(reg, hamiltonian8, 2*pi/10000, π/8, 1)
    @test final_state ≈ statevec(reg1)
end

@testset "operator_evolve_1" begin
    evo_1, final_state = operator_evolve(2π, 1.0, 12.0, 10.0, 1.3, 2.6, 0.35, 10000)
    reg = zero_state(4)
    ω = 1.0
    tau = 2π/ω
    hamiltonian8 = pulse_hamiltonian(ω, 12.0ω , 10.0ω , 1.3ω , 2.6ω , 0.35ω)
    reg1, res1, reg2, res2, evolution_operator = evolve(reg, hamiltonian8, 2*pi/10000, π/8, 1)
    XZZX = XZZX_operator(tau, π/8/tau, 10000)
    evo_block = matblock(evo_1) 
    XZZX_block = matblock(XZZX)
    ff = operator_fidelity(evo_block, XZZX_block)
    @test ff ≈ 0.97
end