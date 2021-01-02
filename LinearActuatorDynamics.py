from sympy import symbols
from sympy.physics.mechanics import *
from sympy.utilities.codegen import make_routine
from sympy.utilities.lambdify import lambdify


## Generating the equations of motion of a double-pendulum system
q1 = dynamicsymbols('q1')
q1d = dynamicsymbols('q1', 1)
u1 = dynamicsymbols('u1')
u1d = dynamicsymbols('u1', 1)
l, m, g = symbols('l m g')

N = ReferenceFrame('N')

O = Point('O')
P = O.locatenew('P', l*N.x)
O.set_vel(N, 0)
P.set_vel(N, u1*N.x)

ParP = Particle('ParP', P, m)

kd = [q1d - u1]
FL = [(P, m*g*N.x)]
BL = [ParP]

KM = KanesMethod(N, q_ind=[q1], u_ind=[u1], kd_eqs=kd)

(fr, frstar) = KM.kanes_equations(BL, FL)
kdd = KM.kindiffdict()
mm = KM.mass_matrix_full
fo = KM.forcing_full
qudots = mm.inv()*fo
qudots = qudots.subs(kdd)
qudots.simplify()
mechanics_printing()
mprint(qudots)

routine = make_routine('derivs', qudots)
print([arg.result_var for arg in routine.results])
print([arg.expr for arg in routine.results])
print([arg.name for arg in routine.arguments])

binary_func = lambdify([u1, g], qudots)
print(binary_func(1, 9.81))

