using FloquetEngineeringToricCode
using Yao
using LinearAlgebra
#定义x，y方向的spin数量及拓扑结构
nx = 2
ny = 2
nq = nx * ny
#function calculate(nq::Int)
#trotterization去逼近含时的调制波形，定义电路layers
nlayers = 100
#定义拓扑结构
topology = periodic_square_lattice(nx, ny)
#定义参数
t = 1.0
dt = t/nlayers
amplitude = 1.0

omega = 2*π
phase = 0.0

Omega_1 = 12.0
Omega_4 = 10.0
J = pi/8
#定义初始态｜0000⟩
state = zero_state(nq)
#initial_amplitudes = zeros(Complex{Float64}, 2^nq)
#initial_amplitudes[0b0000 + 1] = 1.0
#state = ArrayReg(initial_amplitudes)
#根据周期性波形获得分立的强度vector，分别应用在每层电路上
final_state = get_final(t, nq, nlayers, Omega_1, Omega_4, amplitude, omega, phase, topology, state)
final = normalize!(final_state)
state1 = zero_state(nq)
final2 = XZZX_final(t, nq, J, state1)
op = put(nq, 2 => Z)
#stt = statevec(final)
#pauli_z = [1.0 0.0; 0.0 -1.0]
#iden = [1.0 0.0; 0.0 1.0]
#op = kron(iden, pauli_z, iden, iden)
#final_expect = stt' * op * stt
#return final_expect
#end
expectation_values = expect(op, final)
expectation_values2 = expect(op, final2)