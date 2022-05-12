#-----------------------------------------------------------------------
# solution to the compressible Blasius equation (boundary value problem)
# 2(ρμf'')' + ff'' = 0
# (ρμh')' + Prfh' + Pr(γ-1)Ma^{2}ρμf''^{2} = 0
# with isothermal bc f(0) = f'(0) = 0, f'(∞) = 1, h(∞) = 1, h = θ(0).
# following Howarth-Dorodnitsyn transformation θ = T/T_{∞}
#-----------------------------------------------------------------------
using Plots
using Polynomials
using FastGaussQuadrature
using NLsolve

default(linewidth=4)
γ = 1.4  # specific heat ratio
Pr = 0.7 # Prandtl number
Ma = 5.0 # Mach number

# viscosity model for now treat as constants

η_max = 10.0
N = 30
x_resid = gausslegendre(N-3)[1]
dx_dη = 2 / η_max

f = zeros(N)
h = zeros(N)

function CompressibleBlasiusResidual(f, h)
    f = ChebyshevT(f)
    h = ChebyshevT(h)
    f_η = derivative(f, 1) * dx_dη
    f_ηη = derivative(f_η, 1) * dx_dη
    f_ηηη = derivative(f_ηη, 1) * dx_dη
    
end

