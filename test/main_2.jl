using Main.FloquetEngineeringToricCode
using Yao
#定义x，y方向的spin数量及拓扑结构
nx = 2      
ny = 2
nq = nx * ny
topology = periodic_square_lattice(nx, ny)
amplitude = 1.0
omega = 2*π
phase = 0.0
t = 1.0
Omega_1 = 1.0
Omega_4 = 1.0

J = floquet_drive(t, amplitude, omega, phase)


    