using Test, FloquetEngineeringToricCode
using LinearAlgebra
using PauliPropagation
using Yao

@testset "Hamiltonian Evolve" begin
    nq = 2
    q1 = 1
    q2 = 2
    @test mat(xx_yy_time_evolve(nq, q1, q2, 0.0)) ≈ Matrix(LinearAlgebra.I, 4, 4)
    @test mat(xx_yy_time_evolve(nq, q1, q2, 2*π)) ≈ Matrix(LinearAlgebra.I, 4, 4)
    @test isapprox(mat(xx_yy_time_evolve(nq, q1, q2, π / 3)) , 
    [1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
        0.0+0.0im -0.5-0.0im 0.0-0.866025im 0.0+0.0im;
        0.0+0.0im 0.0-0.866025im -0.5-0.0im 0.0+0.0im;
        0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im], atol=1e-6)
end

@testset "trotterization_xxyy evolve" begin
    q1 = 1      
    q2 = 2
    circuit = trotterization_xxyy(q1, q2, PauliString(2, :Z, 2))
    

    
    @test length(circuit) == 400
    
end
@testset "trotterization_circuit" begin
    nx = 2      
    ny = 2
    nq = nx * ny
    
    step = 100
    nlayers = step
    topology = periodic_square_lattice(nx, ny)
    circuit = trotterization_circuit(nq, nlayers; topology=topology)
    # dt: step size
    total_time = 1 
    dt =  total_time / step
    amplitude = 1.0
    omega = 2*π
    phase = 0.0
    Omega_1 = 1.0
    Omega_4 = 1.0
    my_drive = [floquet_drive(cur_t,amplitude,omega,phase) for cur_t in 0:dt:total_time-dt] 
    my_drive_2 = [floquet_drive(cur_t,amplitude,2*omega,phase) for cur_t in 0:dt:total_time-dt]
    my_drive = vcat([[i,i,my_drive_2[idx],my_drive_2[idx],my_drive_2[idx],my_drive_2[idx],Omega_1,Omega_4] for (idx, i) in enumerate(my_drive)]...)
    # my_drive = floquet_drive.(0:dt:total_time-dt,Ref(amplitude),Ref(omega),Ref(phase))

    observable = PauliString(nq, :Z, 2)
    pauli_sum = propagate(circuit,observable,my_drive)
    @test length(circuit) == 800
end

@testset "Circuitwithyao" begin
    nx = 2      
    ny = 2
    nq = nx * ny
    step = 100
    nlayers = step
    t = 1.0
    Omega_1 = 1.0
    Omega_4 = 1.0
    dt = t/nlayers
    topology = periodic_square_lattice(nx, ny)
    coupling_strength = ones(length(topology))
    circuit = circuitwithyao(dt, nq, Omega_1, Omega_4, coupling_strength, topology)
    expected_length =  (length(topology) + 2)
    @test length(circuit) == expected_length
end

function XZZX_circuit(nq::Int, J::Float64, t::Float64)
    state = zero_state(nq)
    circuit = Yao.AbstractBlock[]
    H_coupling = J * kron(X, Z, Z, X)  # X_1 Z_2 Z_3 X_4
    push!(circuit, put(nq, (1,2,3,4) => TimeEvolution(matblock(H_coupling), t)))
    #apply!(state, circuit)
    return circuit
end

  #  state = zero_state(nq)
   # reg = reg |> put(1=>X) |> put(2=>Z) |> put(3=>Z) |> put(4=>X)

function XZZX_final(t::Float64, nq::Int, J::Float64, state::ArrayReg)
    finalstate = zero_state(nq)
    circuit_blocks = XZZX_circuit(nq, J, t)
    circuit = chain(circuit_blocks...)
    finalstate = apply!(state, circuit)
    finalstate = normalize!(finalstate)
    return finalstate
end

