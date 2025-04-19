using Yao
using CairoMakie
using LinearAlgebra
using CairoMakie: PolarAxis

function plot_results(res1, res2; pulses=1:4)
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1], 
        xlabel = "Time", 
        ylabel = "⟨Z₁⟩",
        title = "Comparison of Evolution Methods")
    # Plot both results
    for i in pulses
        lines!(ax, getindex.(res1, i), label = "Pulse Hamiltonian $i", linewidth = 2, color = [:blue, :green, :purple, :orange][i])
        lines!(ax, getindex.(res2, i), label = "XZZX Evolution $i", linewidth = 2, color = [:blue, :green, :purple, :orange][i], linestyle = :dash)
    end
    # Add legend
    axislegend(ax, position = :lb)
    return fig
end

function pulse_hamiltonian(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)
    h(t) = Ω₁ * put(4, 1=>X) + Ω₄ * put(4, 4=>X) +
        g13 * cos(ω * t) * put(4, (1,3)=>xy) +
        g23 * cos(2ω * t) * put(4, (2,3)=>xy) +
        g24 * cos(2ω * t) * put(4, (2,4)=>xy)
    return h
end
function h1(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)    
    h(t) = g13 * cos(ω * t) * put(4, (1,3)=>xy)
    return h
end

function h2(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)    
    h(t) = 
        g23 * cos(2ω * t) * put(4, (2,3)=>xy) 
    return h
end
function h3(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)    
    h(t) = 
    g24 * cos(2ω * t) * put(4, (2,4)=>xy)
    return h
end
function hh(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)    
    h(t) = Ω₁ * put(4, 1=>X) + Ω₄ * put(4, 4=>X) 
    return h
end
function pulse_h(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)    
    h(t) = Ω₁ * put(4, 1=>X) + Ω₄ * put(4, 4=>X) +
        g13 * cos(ω * t) * put(4, (1,2)=>xy) +
        g23 * cos(2ω * t) * put(4, (2,3)=>xy) +
        g24 * cos(2ω * t) * put(4, (3,4)=>xy)
    return h
end

function pulse_hamiltonian_interaction(ω, Ω₁, Ω₄, g13, g23, g24)
    xx = kron(X, X)
    yy = kron(Y, Y)
    zy = kron(Z, Y)
    yz = kron(Y, Z)
    @show Ω₁, Ω₄
    h(t) = g13 * cos(ω * t) * put(4, (1,3)=>xx) +
        g13 * cos(ω * t) * cos(2Ω₁ * t) * put(4, (1,3)=>yy) -
        g13 * cos(ω * t) * sin(2Ω₁ * t) * put(4, (1,3)=>zy) +
        g23 * cos(2ω * t) * put(4, (2,3)=>xx) +
        g23 * cos(2ω * t) * put(4, (2,3)=>yy) +
        g24 * cos(2ω * t) * put(4, (2,4)=>xx) +
        g24 * cos(2ω * t) * cos(2Ω₄ * t) * put(4, (2,4)=>yy) -
        g24 * cos(2ω * t) * sin(2Ω₄ * t) * put(4, (2,4)=>yz)
    return h
end

function pulse_hamiltonian_interaction_1(ω, g13, g23, g24)
    xx = kron(X, X)
    yy = kron(Y, Y)
    
    h(t) = g13 * cos(ω * t) * put(4, (1,3)=>xx) +
        g23 * cos(2ω * t) * put(4, (2,3)=>xx) +
        g23 * cos(2ω * t) * put(4, (2,3)=>yy) +
        g24 * cos(2ω * t) * put(4, (2,4)=>xx) 
    return h
end

function evolve(reg, hamiltonian, t, J, Nt=100000)
    reg1 = copy(reg)
    reg2 = copy(reg)
    dt = t / Nt
    res1 = Vector{Float64}[]
    res2 = Vector{Float64}[]
    hh = -J * (put(4, (1,2,3,4)=>kron(X,Z,Z,X)) + put(4, (1,2)=>kron(Z,X)) + put(4, (3,4)=>kron(X,Z)))
    
    # 计算一个完整周期后的演化算符矩阵
    evolution_operator = Matrix(I, 16, 16)  # 初始化为单位矩阵 
    for it = 1:Nt
        h = hamiltonian((it-0.5) * dt)
        apply!(reg1, time_evolve(h, dt))
        apply!(reg2, time_evolve(hh, dt))
        push!(res1, [expect(put(4, i=>Z), reg1) for i = 1:4])
        push!(res2, [expect(put(4, i=>Z), reg2) for i = 1:4]) 
        evolution_operator = mat(time_evolve(h, dt)) * evolution_operator
    end
    return reg1, res1, reg2, res2, evolution_operator
end

function operator_evolve1(t::Float64, ω::Float64, Ω₁::Float64, Ω₄::Float64, g13::Float64, g23::Float64, g24::Float64, Nt::Int=10000)
    dt = t / Nt
    h_1 = h1(ω, Ω₁, Ω₄, g13, g23, g24)
    h_2 = h2(ω, Ω₁, Ω₄, g13, g23, g24)
    h_3 = h3(ω, Ω₁, Ω₄, g13, g23, g24)
    h_h = hh(ω, Ω₁, Ω₄, g13, g23, g24)
    evolution_operator_k = Matrix(I, 16, 16)  # 初始化为单位矩阵 
    evolution_operator = Matrix(I, 16, 16)  # 初始化为单位矩阵 
    s1 = 1/(4-4^(1/3))
    s2 = (1-4^(1/3))/(4-4^(1/3))
    s3 = 1-2*(s1+s2)
    s4 = s2
    s5 = s1
    s = [s1, s2, s3, s4, s5]
    tt = zeros(5)
    for it = 0:Nt-1
        for k = 1:5
            tt[k] = it*dt + s[k]*dt
            evolution_operator_1 = mat(time_evolve(h_1(tt[k]), s[k]*dt/2)) 
            evolution_operator_2 = mat(time_evolve(h_2(tt[k]), s[k]*dt/2))
            evolution_operator_3 = mat(time_evolve(h_3(tt[k]), s[k]*dt))
            evolution_operator_4 = mat(time_evolve(h_h(tt[k]), s[k]*dt/2))

            evolution_operator_k = (evolution_operator_4 * evolution_operator_1 * evolution_operator_2 * evolution_operator_3 * evolution_operator_2 * evolution_operator_1 * evolution_operator_4) * evolution_operator_k
        end
        evolution_operator = evolution_operator_k * evolution_operator
    end
    return evolution_operator 
end

function operator_evolve(t::Float64, ω::Float64, Ω₁::Float64, Ω₄::Float64, g13::Float64, g23::Float64, g24::Float64, Nt::Int=10000)
    dt = t / Nt
    reg1 = zero_state(4)
    #res1 = Vector{Float64}[]
    h_1 = h1(ω, Ω₁, Ω₄, g13, g23, g24)
    h_2 = h2(ω, Ω₁, Ω₄, g13, g23, g24)
    h_3 = h3(ω, Ω₁, Ω₄, g13, g23, g24)
    h_h = hh(ω, Ω₁, Ω₄, g13, g23, g24)
    evolution_operator = Matrix(I, 16, 16)  # 初始化为单位矩阵 
    for it = 1:Nt
        evolution_operator_1 = mat(time_evolve(h_1((it-0.5) * dt), dt/2)) 
        evolution_operator_2 = mat(time_evolve(h_2((it-0.5) * dt), dt/2))
        evolution_operator_3 = mat(time_evolve(h_3((it-0.5) * dt), dt))
        evolution_operator_4 = mat(time_evolve(h_h((it-0.5) * dt), dt/2))
        evolution_operator = (evolution_operator_4 * evolution_operator_1 * evolution_operator_2 * evolution_operator_3 * evolution_operator_2 * evolution_operator_1 * evolution_operator_4) * evolution_operator


    end
    final_state = evolution_operator * statevec(reg1)
    return evolution_operator, final_state
end
function XZZX_operator(t, J, Nt=10000)
    hh = -J * (put(4, (1,2,3,4)=>kron(X,Z,Z,X)))
    evolution_operator = Matrix(I, 16, 16)  # 初始化为单位矩阵 
    dt = t / Nt
    for it = 1:Nt
        evolution_operator = mat(time_evolve(hh, dt)) * evolution_operator
    end
    return evolution_operator
end



function sim1(ω; interaction=true)
    reg = zero_state(4)
    τ = 2π/ω
    if interaction
        hamiltonian8 = pulse_hamiltonian_interaction(ω, 12.0ω , 10.0ω , 1.3ω , 2.6ω , 0.35ω)
    else
        hamiltonian8 = pulse_hamiltonian(ω, 12.0ω , 10.0ω , 1.3ω , 2.6ω , 0.35ω)
    end
    reg1, res1, reg2, res2, evolution_operator = evolve(reg, hamiltonian8, τ, π/8/τ, 10000)
    return reg1, reg2, res1, res2, evolution_operator
end

function sim2(ω; interaction=false)
    reg = zero_state(4)
    τ = 2π/ω
    if interaction
        hamiltonian50 = pulse_hamiltonian_interaction(ω, 14.0ω, 11.0ω, 1.2ω, 2.6ω, 0.054ω)
    else
        hamiltonian50 = pulse_hamiltonian(ω, 14.0ω, 11.0ω, 1.2ω, 2.6ω, 0.054ω)
    end
    reg1, res1, reg2, res2, evolution_operator = evolve(reg, hamiltonian50, τ, π/50/τ, 20000)
    return reg1, res1, res2, evolution_operator
end

 
# 绘制 Pauli 矩阵迹的极坐标图
function plot_pauli_traces(evo; title="Components of the effective Hamiltonian")
    omega = 1
    tau = 2π/omega

    # 定义 Pauli 矩阵
    XZZX = put(4, (1,2,3,4)=>kron(X,Z,Z,X))
    #YZZY = put(4, (1,2,3,4)=>kron(Y,Z,Z,Y))
    IXIX = put(4, (2,4)=>kron(X, X))
    IYIY = put(4, (2,4)=>kron(Y, Y))
    IZXZ = put(4, (2,3,4)=>kron(Z, X, Z))
    ZZZZ = put(4, (1,2,3,4)=>kron(Z,Z,Z,Z))
    XIXI = put(4, (1,3)=>kron(X, X))
    YIYI = put(4, (1,3)=>kron(Y, Y))
    IXXI = put(4, (2,3)=>kron(X, X))
    IYYI = put(4, (2,3)=>kron(Y, Y))
    ZXZI = put(4, (1,2,3)=>kron(Z,X,Z))
    XIII = put(4, (1)=>X)
    IIIX = put(4, (4)=>X)

    # 计算迹
    traces = [
        tr(evo * (mat(XZZX)))/(16*tau),
        #tr(evo * mat(YZZY))/16,
        tr(evo * (mat(IXIX)))/(16*tau),
        tr(evo * (mat(IYIY)))/(16*tau),
        tr(evo * (mat(IZXZ)))/(16*tau),
        tr(evo * (mat(ZZZZ)))/(16*tau),
        tr(evo * (mat(XIXI)))/(16*tau),
        tr(evo * (mat(YIYI)))/(16*tau),
        tr(evo * (mat(IXXI)))/(16*tau),
        tr(evo * (mat(IYYI)))/(16*tau),
        tr(evo * (mat(ZXZI)))/(16*tau),
        tr(evo * (mat(XIII)))/(16*tau),
        tr(evo * (mat(IIIX)))/(16*tau)
    ]
    @show traces/omega

    # Pauli 矩阵标签
    labels = ["XZZX", "IXIX", "IYIY", "IZXZ", "ZZZZ", "XIXI", 
              "YIYI", "IXXI", "IYYI", "ZXZI", "XIII", "IIIX"]

    # 创建图形
    fig = Figure()
    ax = PolarAxis(fig[1, 1], 
        title = "Pauli Matrix Traces (First Simulation)",
        thetaticks = 0:π/6:2π  # 每30度（π/6弧度）一条径向线
    )
    angles = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]  
    theta = deg2rad.(angles)
    r = log10.(abs.(real.(traces)/omega))
    
    # 绘制散点图
    scatter = scatter!(ax, theta, r, 
        color = :blue,
    )
    
    # 添加标签
    for (i, (angle, label)) in enumerate(zip(theta, labels))
        text!(ax, angle, maximum(r) - 1.0,  # 在最大半径外一点的位置添加标签
            text = label,
            align = (:center, :center),
            fontsize = 13
        )
    end
    
    return fig

end



