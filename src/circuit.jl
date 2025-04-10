#Circuit constructed with Yao.jl
function circuitwithyao(t::Float64, nq::Int, nlayers::Int, Omega_1::Float64, Omega_4::Float64, coupling_strength::Vector{Float64}, topology::Vector{Tuple{Int,Int}})
    dt = t/nlayers
    circuit = chain(nq,
    mapreduce(vcat,1:nlayers) do i 
       [[xx_yy_time_evolve(nq, pair..., coupling_strength[j]*dt/5) for (j, pair) in enumerate(topology)]; [rot(put(nq,2=>X),Omega_1*dt/5), rot(put(nq,3=>X),Omega_4*dt/5)]]
    end 
    )
    return circuit
end

# modularize,make time evolution of xx+yy separately
xx_yy_time_evolve(nq::Int, q1::Int,q2::Int, dt::Float64) = time_evolve(kron(nq, q1=>X, q2=>X) + kron(nq, q1=>Y, q2=>Y), dt) 


#trotterization circuit with n layers
function trotterization_circuit(nq::Int, nlayers::Int; topology::Vector{Tuple{Int,Int}})
    circuit::Vector{Gate} = []
    for i in 1:nlayers
        for cur_topo in topology
            push!(circuit, PauliRotation([:X, :X], cur_topo))
            push!(circuit, PauliRotation([:Y, :Y], cur_topo))
        end
        push!(circuit, PauliRotation(:X, 2))
        push!(circuit, PauliRotation(:X, 3))
    end
    return circuit
end

# evolution of XX+YY separately
function trotterization_xxyy(q1::Int, q2::Int, Observable::PauliString)
    circuit::Vector{Gate} = []
    push!(circuit, PauliRotation([:X, :X], (q1,q2)))
    push!(circuit, PauliRotation([:Y, :Y], (q1,q2)))
    pauli_sum = propagate(circuit, Observable, [1.0,1.0])
    return pauli_sum
end
