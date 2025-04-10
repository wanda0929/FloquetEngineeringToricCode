using Test, FloquetEngineeringToricCode
using PauliPropagation
using Yao

@testset "periodic_square_lattice" begin
    # smallest size of the lattice that has non trivial lattice
    # for all three shifts
    # lattice indices
    # 1 2 3
    # 4 5 6
    # 7 8 9
    nx = 3
    ny = 3
    nq = nx * ny
    lattice_00 = periodic_square_lattice(nx,ny)
    lattice_01 = periodic_square_lattice(nx,ny,xshift=1)
    lattice_10 = periodic_square_lattice(nx,ny,yshift=1)
    lattice_11 = periodic_square_lattice(nx,ny,xshift=1,yshift=1)
    @test lattice_00 ==  [(1,4),(1,5),(2,5)]
    @test lattice_01 ==  [(2,5),(2,6),(3,6)]
    @test lattice_10 ==  [(4,7),(4,8),(5,8)]
    @test lattice_11 ==  [(5,8),(5,9),(6,9)]
end


@testset "circuitofpaulipropogation" begin
    nx = 2      
    ny = 2
    nq = nx * ny
    step = 1
    nlayers = step
    dt = 1.0/nlayers
    topology = periodic_square_lattice(nx, ny)
    circuit = trotterization_circuit(nq, nlayers; topology=topology)
    circuit_1 = circuitwithyao(dt, nq, 10.0, 10.0, [1.0, 1.0, 1.0], topology)
    println(circuit)
    state = zero_state(nq)
    observable = PauliString(nq, :Z, 2)
    pauli_sum = propagate(circuit,observable,[1.0,1.0,1.0,1.0,1.0,1.0,10.0,10.0])
    final_pauli = overlapwithzero(pauli_sum)
    circuit = chain(circuit_1...)
    finalstate = apply!(state, circuit)
    final = normalize!(finalstate)
    op = put(nq, 2 => Z)
    expectation_values = expect(op, finalstate)
    @test final_pauli == expectation_values
end


@testset "Getfinal" begin
    nx = 2      
    ny = 2
    nq = nx * ny
    step = 100
    nlayers = step
    t = 1.0
    Omega_1 = 1.0
    Omega_4 = 1.0
    topology = periodic_square_lattice(nx, ny)
    initial_amplitudes = zeros(Complex{Float64}, 2^nq)
    initial_amplitudes[0b0000 + 1] = 1.0
    state = ArrayReg(initial_amplitudes)
    get_final(t, nq, nlayers, Omega_1, Omega_4, 1.0, 2*π, 0.0, topology, state)
end

@testset "XZZX_circuit" begin
    nx = 2      
    ny = 2
    nq = nx * ny
    J = 1.0
    t = 0.1
    circuit = XZZX_circuit(nq, J, dt)
    @test length(circuit) == 1
end

@testset "XZZX_final" begin
    nx = 2      
    ny = 2
    nq = nx * ny
    J = 1.0
    t = 0.1
    circuit = XZZX_circuit(nq, J, dt)
    initialstate = zero_state(nq)
    final = XZZX_final(t, nq, J, initialstate)


end