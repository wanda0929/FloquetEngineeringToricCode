#定义系统的拓扑结构
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

#定义用Yao.jl构建的电路
#function circuitwithyao(dt::Float64, nq::Int, Omega_1::Float64, Omega_4::Float64, coupling_strength::Vector{Float64}, topology::Vector{Tuple{Int,Int}})
#    circuit = Yao.AbstractBlock[]
    #dt = t/nlayers
    #for i in 1:nlayers
       
#        for (j, pair) in enumerate(topology)
#           J = coupling_strength[j]
#            H_1 = J * kron(X, X)
#            H_2 = J * kron(Y, Y)
#            push!(circuit, put(nq, pair => TimeEvolution(matblock(H_1), dt)))
#            push!(circuit, put(nq, pair => TimeEvolution(matblock(H_2), dt)))
#        end
#        push!(circuit, put(nq, 2 => TimeEvolution(matblock(Omega_1 * X), dt)))
#        push!(circuit, put(nq, 3 => TimeEvolution(matblock(Omega_4 * X), dt)))
    #end
#    return circuit
#end


#定义Yao.jl分层应用在初始态上的函数
function get_final(t::Float64, nq::Int, nlayers::Int, Omega_1::Float64, Omega_4::Float64, amplitude::Float64, omega::Float64, phase::Float64, topology::Vector{Tuple{Int,Int}},state::ArrayReg)
    dt = t/nlayers
    finalstate = zero_state(nq)
    for i in 1:nlayers
        time_t = i*dt
        J_1 = floquet_drive(time_t, amplitude, omega, phase)
        J_2 = floquet_drive(time_t, amplitude, 2*omega, phase)
        coupling_strength = [J_1 ,J_2, J_2]
        circuit_blocks = circuitwithyao(dt, nq, Omega_1, Omega_4, coupling_strength, topology)
        circuit = chain(circuit_blocks...)
        finalstate = apply!(state, circuit)
    end
    return finalstate
end
"""
function initial_state(nq::Int)
    circuit = Yao.AbstractBlock[]
    iden = [1.0 0.0; 0.0 1.0]
    H = 0.5*(kron((kron(iden,iden)+kron(Z,X)),(kron(iden,iden) + kron(X,Z))))
    op_block = matblock(H)
    state = zero_state(nq)
    reg = reg |> put(1=>H) |> put(1=>Z) |> put(2=>X) |> put(3=>H) |> put(3=>Z) |> put(4=>X)

    initial_amplitudes = zeros(Complex{Float64}, 2^nq)
    initial_amplitudes[0b0000 + 1] = 1.0
    return ArrayReg(initial_amplitudes)
end
"""
