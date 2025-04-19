using Test
using CairoMakie
using LinearAlgebra
using FloquetEngineeringToricCode
using Yao, Yao.EasyBuild
using KrylovKit: eigsolve
using Optimisers
using LinearAlgebra
#using DifferentialEquations

@testset "pulse_hamiltonian_seperated" begin
    t = 1.0
    ω = 1
    Ω₁ = 1.0
    Ω₄ = 1.0
    params = [1.0, 1.0, 1.0]
    ham1, ham2, ham3, ham4 = pulse_hamiltonian_seperated(ω, Ω₁, Ω₄, params)
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
    Ω₁ = 1.0
    Ω₄ = 1.0
    params = [1.0, 1.0, 1.0]
    operator = operator_evolve(t, ω, Ω₁, Ω₄, params)
    final_state = operator * initial
    pulse_ham = pulse_hamiltonian0(ω, Ω₁, Ω₄, params[1], params[2], params[3])
    reg1, res1, reg2, res2, evolution_operator = evolve(initial_state, pulse_ham, t, 1.0)
    @test statevec(reg1) ≈ final_state
end

@testset "XZZX_operator" begin
    t = 1.0
    J = 1.0
    operator = XZZX_operator(t, J)
    @test operator isa Yao.AbstractMatrix
end

@testset "diff_evolution_operator" begin
    t = 1.0
    J = 1.0
    ω = 1
    Ω₁ = 1.0
    Ω₄ = 1.0
    params = [1.0, 1.0, 1.0]
    operator1 = XZZX_operator(t, J)
    operator2 = operator_evolve(t, ω, Ω₁, Ω₄, params)
    diff = operator1 - operator2
    diff_evolution_operator(t, ω, Ω₁, Ω₄, params, J)
    @test norm(diff) ≈ diff_evolution_operator(t, ω, Ω₁, Ω₄, params, J)
end

@testset "gradient" begin
    t = 1.0
    J = 1.0
    ω = 1
    Ω₁ = 1.0
    Ω₄ = 1.0
    params = [1.0, 1.0, 1.0]
    grad = gradient(t, ω, Ω₁, Ω₄, params, J)
    @test grad isa Vector{Float64}
end

@testset "optimize_process" begin
    t = 1.0
    J = 1.0
    ω = 1
    Ω₁ = 0.0
    Ω₄ = 0.0
    params = [1.0, 1.0, 1.0]
    Nt = 10000
    ε = 1e-6
    optimized = optimize_process(t, ω, Ω₁, Ω₄, J)
    @test optimized isa Vector{Float64}
end


