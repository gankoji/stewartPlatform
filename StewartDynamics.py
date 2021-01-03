from sympy import symbols
from sympy.physics.vector.point import Point
from sympy.physics.mechanics import *
from sympy.utilities.codegen import make_routine
from sympy.utilities.lambdify import lambdify
from sympy.matrices import Matrix, eye

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

## Coordinates of actuators on moving platform
x1, y1 = symbols('x1 y1')
x2, y2 = symbols('x2 y2')
x3, y3 = symbols('x3 y3')
x4, y4 = symbols('x4 y4')
x5, y5 = symbols('x5 y5')
x6, y6 = symbols('x6 y6')

## Coordinates of actuators on base platform
bx1, by1 = symbols('bx1 by1')
bx2, by2 = symbols('bx2 by2')
bx3, by3 = symbols('bx3 by3')
bx4, by4 = symbols('bx4 by4')
bx5, by5 = symbols('bx5 by5')
bx6, by6 = symbols('bx6 by6')

## Geometry specification
N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)

## Attachment points on base platform
B1 = O.locatenew('B1', x1*N.x + y1*N.y)
B2 = O.locatenew('B2', x2*N.x + y2*N.y)
B3 = O.locatenew('B3', x3*N.x + y3*N.y)
B4 = O.locatenew('B4', x4*N.x + y4*N.y)
B5 = O.locatenew('B5', x5*N.x + y5*N.y)
B6 = O.locatenew('B6', x6*N.x + y6*N.y)


## Platform body frame and inertia
B = ReferenceFrame('B')
B.orient(N, 'DCM', eye(3), 'ZYX')
I = inertia(B, Ixx, Iyy, Izz)

## Center of mass of moving platform
P = O.locatenew('P', h*N.z)

## Attachment points on moving platform
T1 = P.locatenew('T1', x1*B.x + y1*B.y)
T2 = P.locatenew('T2', x2*B.x + y2*B.y)
T3 = P.locatenew('T3', x3*B.x + y3*B.y)
T4 = P.locatenew('T4', x4*B.x + y4*B.y)
T5 = P.locatenew('T5', x5*B.x + y5*B.y)
T6 = P.locatenew('T6', x6*B.x + y6*B.y)

Plat = RigidBody('Plat', P, B, m, (I,P))

const_coord = [T1.pos_from(B1).magnitude() - q1, T2.pos_from(B2).magnitude() - q2, T3.pos_from(B3).magnitude() - q3, T4.pos_from(B4).magnitude() - q4, T5.pos_from(B5).magnitude() - q5, T6.pos_from(B6).magnitude() - q6]

kd = [q1d - u1, q2d - u2, q3d - u3, q4d - u4, q5d - u5, q6d - u6]
FL = [(P, -m*g*N.z)]
BL = [Plat]

KM = KanesMethod(N, q_ind=[q1, q2, q3, q4, q5, q6], u_ind=[u1, u2, u3, u4, u5, u6], kd_eqs=kd)
               ## configuration_constraints=const_coord)

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

