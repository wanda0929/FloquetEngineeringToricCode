using Yao, Yao.EasyBuild
using KrylovKit
import KrylovKit: eigsolve
using Optimisers
using CairoMakie
using LinearAlgebra
import LinearAlgebra: norm

function pulse_hamiltonian0(ω, Ω₁, Ω₄, g13, g23, g24)
    xy = kron(X, X) + kron(Y, Y)
    h(t) = Ω₁ * put(4, 1=>X) + Ω₄ * put(4, 4=>X) +
        g13 * cos(ω * t) * put(4, (1,3)=>xy) +
        g23 * cos(3ω * t) * put(4, (2,3)=>xy) +
        g24 * cos(2ω * t) * put(4, (2,4)=>xy)
    return h
end
function pulse_hamiltonian_seperated(ω, Ω₁, Ω₄, params)
    xy = kron(X, X) + kron(Y, Y)
    coeff1(t) = params[1] * cos(ω * t)
    coeff2(t) = params[2] * cos(3ω * t)
    coeff3(t) = params[3] * cos(2ω * t)

    # 预计算基础矩阵
    base_ham1 = mat(Ω₁ * put(4, 1=>X) + Ω₄ * put(4, 4=>X))
    base_ham2 = mat(put(4, (1,3)=>xy))
    base_ham3 = mat(put(4, (2,3)=>xy))
    base_ham4 = mat(put(4, (2,4)=>xy))

    # 返回矩阵函数
    ham1 = base_ham1
    ham2(t) = coeff1(t) * base_ham2
    ham3(t) = coeff2(t) * base_ham3
    ham4(t) = coeff3(t) * base_ham4
    return ham1, ham2, ham3, ham4
end

function operator_evolve(t, ω, Ω₁, Ω₄, params, Nt=10000)
    dt = t / Nt
    ham1, ham2, ham3, ham4 = pulse_hamiltonian_seperated(ω, Ω₁, Ω₄, params)
    evolution_operator = Matrix{ComplexF64}(I, 16, 16)  # 初始化为单位矩阵 
    for it = 1:Nt
        t_current = (it-0.5) * dt
        evolution_operator_1 = exp(-im * (dt/2) * Matrix(ham1))
        evolution_operator_2 = exp(-im * (dt/2) * Matrix(ham2(t_current)))
        evolution_operator_3 = exp(-im * (dt/2) * Matrix(ham3(t_current)))
        evolution_operator_4 = exp(-im * (dt) * Matrix(ham4(t_current)))
        evolution_operator = (evolution_operator_1 * (evolution_operator_2 * evolution_operator_3 * evolution_operator_4 * evolution_operator_3 * evolution_operator_2) * evolution_operator_1) * evolution_operator
    end
    return evolution_operator
end

function XZZX_operator(t, J)
    hh = -J * (put(4, (1,2,3,4)=>kron(X,Z,Z,X))) |> mat
    evolution_operator = exp(-im * t * hh)
    return evolution_operator
end

function diff_evolution_operator(t, ω, Ω₁, Ω₄, params, J, Nt=10000)
    op_1 = operator_evolve(t, ω, Ω₁, Ω₄, params, Nt)
    op_2 = XZZX_operator(t, J)
    norm_op = norm(op_1 - op_2)
    return norm_op
end

function gradient(t, ω, Ω₁, Ω₄, params, J, Nt=10000, ε=1e-6)
    grad = zeros(length(params))
    diff_0 = diff_evolution_operator(t, ω, Ω₁, Ω₄, params, J, Nt)
    for i in 1:length(params)
        params_new = copy(params)
        params_new[i] += ε
        diff_1 = diff_evolution_operator(t, ω, Ω₁, Ω₄, params_new, J, Nt)
        grad[i] = (diff_1 - diff_0) / ε
    end
    return grad
end


function optimize_process(t, ω, Ω₁, Ω₄, J, Nt=10000, ε=1e-6)
    params = rand(3)
    optimizer = Optimisers.setup(Optimisers.ADAM(0.01), params)
    niter = 100
    for i = 1:niter
        grad_params = gradient(t, ω, Ω₁, Ω₄, params, J, Nt, ε)
        optimizer, params = Optimisers.update(optimizer, params, grad_params)
        current_loss = diff_evolution_operator(t, ω, Ω₁, Ω₄, params, J, Nt)
        println("Iteration $i: Loss = $current_loss")
    end
    return params
end

