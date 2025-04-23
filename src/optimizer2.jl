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
        g23 * cos(2ω * t) * put(4, (2,3)=>xy) +
        g24 * cos(2ω * t) * put(4, (2,4)=>xy)
    return h
end
function pulse_hamiltonian_seperated(ω, params,param_optimized)
    xy = kron(X, X) + kron(Y, Y)
    Ω₁ = params[1]
    Ω₄ = params[2]
    coeff1(t) = (param_optimized[1] + 0.01 * params[3]) * cos(1 * ω * t)
    coeff2(t) = (param_optimized[2] + 0.01 * params[4]) * cos(2 * ω * t)
    coeff3(t) = (param_optimized[3] + 0.01 * params[5]) * cos(2 * ω * t)

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

function operator_evolve(t, ω, params, param_optimized, Nt=10000)
    dt = t / Nt
    ham1, ham2, ham3, ham4 = pulse_hamiltonian_seperated(ω, params,param_optimized)
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

function diff_evolution_operator(t, ω, params, param_optimized, J, Nt=10000)
    op_1 = operator_evolve(t, ω, params, param_optimized, Nt)
    op_2 = XZZX_operator(t, J)
    norm_op = norm(op_1 - op_2)
    return norm_op
end

function gradient(t, ω, params, param_optimized, J, Nt=10000, ε=1e-3)
    grad = zeros(length(params))
    diff_0 = diff_evolution_operator(t, ω, params, param_optimized, J, Nt)
    for i in 1:length(params)
        params_new = copy(params)
        params_new[i] += ε
        diff_1 = diff_evolution_operator(t, ω, params_new, param_optimized, J, Nt)
        grad[i] = (diff_1 - diff_0) / ε
    end
    return grad
end

# 参数变换函数
function transform_params(log_params)
    return exp.(log_params)  # 将参数从对数空间转换回原始空间
end

function inverse_transform_params(params)
    return log.(params)  # 将参数转换到对数空间
end

function optimize_process(t, ω, J, initial_params, param_optimized, Nt=10000, ε=1e-3)
    if any(x -> x < 0, initial_params)
        error("All initial parameters must be positive")
    end
    log_params = inverse_transform_params(initial_params)
    println("Initial log parameters: ", log_params)
    optimizer = Optimisers.setup(Optimisers.ADAM(0.01), log_params)
    niter = 100

    best_loss = Inf
    best_params = nothing

    for i = 1:niter
        # 计算梯度
        current_params = transform_params(log_params)
        grad_log_params = gradient(t, ω, current_params, param_optimized, J, Nt, ε)
        
        # 更新参数
        optimizer, log_params = Optimisers.update(optimizer, log_params, grad_log_params)
        
        # 计算当前loss
        current_params = transform_params(log_params)
        current_loss = diff_evolution_operator(t, ω, current_params, param_optimized, J, Nt)
        
        # 更新最佳结果
        if current_loss < best_loss
            best_loss = current_loss
            best_params = copy(current_params)
        end

    println("Best loss achieved: $best_loss")
    println("Best parameters: $best_params")
    end
    return best_params, best_loss
end

function optimize_with_multistarts(t, ω, J, param_optimized, n_starts=10)
    best_loss = Inf
    best_optimized = nothing
    for i = 1:n_starts
        initial_params = 15*rand(5)
        optimized, loss = optimize_process(t, ω, J, initial_params, param_optimized)
        if loss < best_loss
            best_loss = loss
            best_optimized = optimized
        end
    end
    return best_optimized, best_loss
end
