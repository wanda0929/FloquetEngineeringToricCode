function periodic_square_lattice(nx::Int, ny::Int;xshift::Int=0,yshift::Int=0)
    topology = Tuple{Int,Int}[]
    for i in 1+xshift:2:nx-1
        for j in 1+yshift:2:ny-1
            cur_pt = (j-1)*nx + i
            rt_pt = i == nx ? nx*ny + 1  : (j-1)*nx + mod1(i+1,nx)
            dw_pt = j >= ny ? nx*ny + 1 : mod1(j,ny)*nx + i
            rt_dw_pt =  (j >= ny || i == nx ) ? nx*ny+1 : mod1(j,ny)*nx + mod1(i+1,nx)
            dw_pt <= ny*ny && push!(topology, tuple(sort([cur_pt, dw_pt])...))    
            rt_dw_pt <= ny*ny && push!(topology,  tuple(sort([cur_pt, rt_dw_pt])...))
            (rt_dw_pt <= nx*ny && rt_pt <= nx*ny ) && push!(topology, tuple(sort([rt_dw_pt, rt_pt])...))
        end
    end
    return unique(topology)
end

function floquet_drive(t::Float64, amplitude::Float64, omega::Float64, phase=0.0)
    return amplitude * cos(omega * t + phase)
end

 
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

function circuitwithyao(t::Float64, nq::Int, nlayers::Int, Omega_1::Float64, Omega_4::Float64, coupling_strength::Vector{Float64}, topology::Vector{Tuple{Int,Int}})
    circuit = Yao.AbstractBlock[]
    dt = t/nlayers
    for i in 1:nlayers
       
        for (j, pair) in enumerate(topology)
            J = coupling_strength[j]
            H = J * (kron(X, X) + kron(Y, Y))
            push!(circuit, put(nq, pair => TimeEvolution(matblock(H), dt)))
        end
        push!(circuit, put(nq, 2 => TimeEvolution(matblock(Omega_1 * X), dt)))
        push!(circuit, put(nq, 3 => TimeEvolution(matblock(Omega_4 * X), dt)))
    end
    return circuit
end



