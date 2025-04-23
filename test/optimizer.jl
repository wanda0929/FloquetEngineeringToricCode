using Test
using CairoMakie
using LinearAlgebra
using FloquetEngineeringToricCode
using Yao, Yao.EasyBuild
using KrylovKit: eigsolve
using Optimisers
using LinearAlgebra
#using DifferentialEquations

@testset "pulse_hamiltonian_seperated0" begin
    t = 1.0
    ω = 1
    params = [1.0, 1.0, 1.0]
    ham1, ham2, ham3, ham4 = pulse_hamiltonian_seperated0(ω, params)
    ham = pulse_hamiltonian00(ω, params[1], params[2], params[3])
    #@test ham isa Yao.AbstractMatrix
    #@test ham2 isa Yao.AbstractMatrix
    #@test ham3 isa Yao.AbstractMatrix
    #@test ham4 isa Yao.AbstractMatrix
    Matrix(ham(1.0)) ≈ (Matrix(ham1) + Matrix(ham2(1.0)) + Matrix(ham3(1.0)) + Matrix(ham4(1.0)))
end

@testset "operator_evolve0" begin
    initial_state = zero_state(4)
    initial = statevec(initial_state)
    t = 1.0
    ω = 1
    params = [1.0, 1.0, 1.0]
    operator = operator_evolve0(t, ω, params)
    final_state = operator * initial
    pulse_ham = pulse_hamiltonian00(ω, params[1], params[2], params[3])
    reg1, res1, reg2, res2, evolution_operator = evolve(initial_state, pulse_ham, t, 1.0)
    @test statevec(reg1) ≈ final_state
end

@testset "XZZX_operator" begin
    t = 1.0
    J = 1.0
    operator = XZZX_operator(t, J)
    @test operator isa Yao.AbstractMatrix
end

@testset "diff_evolution_operator0" begin
    t = 1.0
    J = 1.0
    ω = 1
    params = [1.0, 1.0, 1.0]
    operator1 = XZZX_operator(t, J)
    operator2 = operator_evolve0(t, ω, params)
    diff = operator1 - operator2
    diff_evolution_operator0(t, ω, params, J)
    @test norm(diff) ≈ diff_evolution_operator0(t, ω, params, J)
end

@testset "gradient0" begin
    t = 1.0
    J = 1.0
    ω = 1
    params = [1.0, 1.0, 1.0]
    grad = gradient0(t, ω, params, J)
    @test grad isa Vector{Float64}
end

@testset "optimize_process0" begin
    J = 1.0
    ω = 1
    t = 2*pi/ω
    params = [1.0, 1.0, 1.0]
    Nt = 10000
    ε = 1e-6
    optimized = optimize_process0(t, ω, J, params)
    @test optimized isa Vector{Float64}
end

@testset "optimize_with_multistarts0" begin
    J = 1.0
    ω = 1
    t = 2*pi/ω
    Nt = 10000
    ε = 1e-6
    optimized, loss = optimize_with_multistarts0(t, ω, J)
    @test optimized isa Vector{Float64}
    @test loss < 10
end




@testset "pulse_hamiltonian_seperated" begin
    t = 1.0
    ω = 1
    Ω₁ = 1.0
    Ω₄ = 1.0
    params = [1.0, 1.0, 1.0, 1.0, 1.0]
    param_optimized = [1.0, 1.0, 1.0]
    ham1, ham2, ham3, ham4 = pulse_hamiltonian_seperated(ω, params,param_optimized)
    ham = pulse_hamiltonian0(ω, Ω₁, Ω₄, params[1], params[2], params[3])
    #@test ham isa Yao.AbstractMatrix
    #@test ham2 isa Yao.AbstractMatrix
    #@test ham3 isa Yao.AbstractMatrix
    #@test ham4 isa Yao.AbstractMatrix
    Matrix(ham(1.0)) ≈ (Matrix(ham1) + Matrix(ham2(1.0)) + Matrix(ham3(1.0)) + Matrix(ham4(1.0)))
end

@testset "operator_evolve" begin
    initial_state = zero_state(4)
    initial = statevec(initial_state)
    t = 1.0
    ω = 1
    params = [1.0, 1.0, 1.0, 1.0, 1.0]
    param_optimized = [1.0, 1.0, 1.0]
    operator = operator_evolve(t, ω, params, param_optimized)
    final_state = operator * initial
    pulse_ham = pulse_hamiltonian0(ω, params[1], params[2], (param_optimized[1] + 0.001 * params[3]), (param_optimized[2] + 0.001 * params[4]), (param_optimized[3] + 0.001 * params[5]))
    reg1, res1, reg2, res2, evolution_operator = evolve(initial_state, pulse_ham, t, 1.0)
    @test statevec(reg1) ≈ final_state
end


@testset "diff_evolution_operator" begin
    t = 1.0
    J = 1.0
    ω = 10
    params = [1.0, 1.0, 1.0, 1.0, 1.0]
    param_optimized = [1.0, 1.0, 1.0]
    operator1 = XZZX_operator(t, J)
    operator2 = operator_evolve(t, ω, params, param_optimized)
    diff = operator1 - operator2
    diff_evolution_operator(t, ω, params, param_optimized, J)
    @test norm(diff) ≈ diff_evolution_operator(t, ω, params, param_optimized, J)
end

@testset "gradient" begin
    t = 1.0
    J = 1.0
    ω = 10
    params = [1.0, 1.0, 1.0, 1.0, 1.0]
    param_optimized = [1.0, 1.0, 1.0]
    grad = gradient(t, ω, params, param_optimized, J)
    @test grad isa Vector{Float64}
end

@testset "optimize_process" begin
    J = 1.0
    ω = 10
    t = 2*pi/ω
    params = [1.0, 1.0, 1.0, 1.0, 1.0]
    param_optimized = [1.0, 1.0, 1.0]
    Nt = 10000
    ε = 1e-6
    optimized = optimize_process(t, ω, J, params, param_optimized)
    @test optimized isa Vector{Float64}
end

@testset "optimize_with_multistarts" begin
    J = 1.0
    ω = 10
    t = 2*pi/ω
    param_optimized = [1.0, 1.0, 1.0]
end

