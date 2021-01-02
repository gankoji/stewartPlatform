from sympy import symbols
from sympy.physics.mechanics import *
from sympy.utilities.codegen import make_routine
from sympy.utilities.lambdify import lambdify


## Generating the equations of motion of a Gough-Stewart Platform system

## Generalized coordinates and speeds
## Here we use the length of the 6 actuators
q1, q2, q3, q4, q5, q6 = dynamicsymbols('q1 q2 q3 q4 q5 q6')
q1d, q2d, q3d, q4d, q5d, q6d = dynamicsymbols('q1 q2 q3 q4 q5 q6', 1)
u1, u2, u3, u4, u5, u6 = dynamicsymbols('u1 u2 u3 u4 u5 u6')
u1d, u2d, u3d, u4d, u5d, u6d = dynamicsymbols('u1 u2 u3 u4 u5 u6', 1)

## Platform parameters
h, m, g = symbols('h m g')
Ixx, Iyy, Izz = symbols('Ixx Iyy Izz')

## Attachment points of actuators on moving platform
x1, y1 = symbols('x1 y1')
x2, y2 = symbols('x2 y2')
x3, y3 = symbols('x3 y3')
x4, y4 = symbols('x4 y4')
x5, y5 = symbols('x5 y5')
x6, y6 = symbols('x6 y6')

## Geometry specification
N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)

B = ReferenceFrame('B')

## Center of mass of moving platform
P = O.locatenew('P', h*N.z)
#P.set_vel(N, u1*N.x)
ParP = Particle('ParP', P, m)

kd = [q1d - u1, q2d - u2, q3d - u3, q4d - u4, q5d - u5, q6d - u6]
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

