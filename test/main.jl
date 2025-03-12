using Main.FloquetEngineeringToricCode
using PauliPropagation
#定义x，y方向的spin数量及拓扑结构
nx = 2      
ny = 2
nq = nx * ny
#trotterization去逼近含时的调制波形，定义电路layers
step = 100
nlayers = step
#定义拓扑结构
topology = periodic_square_lattice(nx, ny)
#定义线路
circuit = trotterization_circuit(nq, nlayers; topology=topology)# dt: step size
#定义参数
total_time = 1 
dt =  total_time / step
amplitude = 1.0
omega = 2*π
phase = 0.0
Omega_1 = 1.0
Omega_4 = 1.0
#根据周期性波形获得分立的强度vector，分别应用在每层电路上
my_drive = [floquet_drive(cur_t,amplitude,omega,phase) for cur_t in 0:dt:total_time-dt] 
my_drive_2 = [floquet_drive(cur_t,amplitude,2*omega,phase) for cur_t in 0:dt:total_time-dt]
my_drive = vcat([[i,i,my_drive_2[idx],my_drive_2[idx],my_drive_2[idx],my_drive_2[idx],Omega_1,Omega_4] for (idx, i) in enumerate(my_drive)]...)
# my_drive = floquet_drive.(0:dt:total_time-dt,Ref(amplitude),Ref(omega),Ref(phase))
#定义observable
#observable = PauliString(nq, :Z, 2)
observable = PauliSum([
    #PauliString(nq, [], []),  # 1
    PauliString(nq, [:Z, :X], [1, 2]),  # Z_1X_2
    PauliString(nq, [:X, :Z], [3, 4]),  # X_3Z_4
    PauliString(nq, [:Z, :X, :X, :Z], [1, 2, 3, 4])  # Z_1X_2X_3Z_4
])
#PauliPropagation
pauli_sum = propagate(circuit,observable,my_drive)