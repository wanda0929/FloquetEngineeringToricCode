
using CairoMakie
using LinearAlgebra
using FloquetEngineeringToricCode
using Yao, Yao.EasyBuild
using KrylovKit: eigsolve
using Optimisers
using LinearAlgebra
using Distributions


ω = 1
t = 2*pi/ω
J = 8/(pi*2*pi/ω)
#params = [0.001, 0.001, 0.14, 1.76, 0.32, 1.77, 3.14, 2.74]
optimized = optimize_with_multistarts0(t, ω, J)

params = 10*rand(5)
optimizedparams = [-0.9334819667898443, -1.1777334305857767, -0.02868486747271791]
op2 = optimize_with_multistarts(t, ω, J, optimizedparams)




[2.649994597208409, 1.5397548366856941, 0.2493044339493439]

([13.338451046817935, 12.48976204959711, 24.60787267028943, 16.256694923972155, 25.061021486754587], 2.7839124389333962)