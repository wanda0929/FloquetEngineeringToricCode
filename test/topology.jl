using Test, FloquetEngineeringToricCode
using PauliPropagation
using Yao
@testset "periodic_square_lattice" begin
    nx = 5
    ny = 5
    nq = nx * ny
    @test length(periodic_square_lattice(nx, ny)) == 12

    
end

@testset "FloquetDrive" begin
    t = 0.5
    amplitude = 3.0
    frequency = 2.0
    phase = 0
    @test floquet_drive(t, amplitude, frequency,phase) == 3.0 * cos(2.0 * 0.5)
end

@testset "Trotterizationcircuit" begin
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
    topology = periodic_square_lattice(nx, ny)
    coupling_strength = ones(length(topology))
    circuit = circuitwithyao(t, nq, nlayers, Omega_1, Omega_4, coupling_strength, topology)
    expected_length = nlayers * (length(topology) + 2)
    @test length(circuit) == expected_length
end 