read "/path/to/SophiaV6.txt":
read "/path/to/Graphics.txt":
with(plots): with(plottools): with(linalg): with(ArrayTools):

# --------------Declaring time depencece---------------
dependsTime(q1, q2, q3, u1, u2, u3):

# --------------Coordinate Systems---------------------
&rot([N, A, 3, q2]):

# --------------Position vector to balls---------------
r_B1 := (A &ev [R*cos(H0-l0/R-q1),0,R*sin(H0-l0/R-q1)]):
r_B2 := (A &ev [R*cos(H1+l1/R+q3),0,R*sin(H1+l1/R+q3)]):

# ---------Velocity vectors of the balls---------------
v_B1 := simplify(N &fdt r_B1):
v_B2 := simplify(N &fdt r_B2):

# ---------------Angular momentum (H)------------------
omega := <0,0,q2t>:
I_hoop := << ((3/4)*r^2+R^2)*m1 | 0 | 0 >,
           < 0 | ((5/8)*r^2+(1/2)*R^2)*m1 | 0 >,
           < 0 | 0 | ((5/8)*r^2+(1/2)*R^2)*m1 >>:
H_hoop := I_hoop.omega:
H_B1 := Ecross(r_B1,m0 &**v_B1):
H_B2 := Ecross(r_B2,m0 &**v_B2):

# ----------------Kinetic energy (T)-------------------
T_hoop_rotation := (1/2)*(omega^+).H_hoop:
T_B1_translation := (1/2)*m0*(v_B1 &o v_B1):
T_B2_translation := (1/2)*m0*(v_B2 &o v_B2):
T := T_hoop_rotation + T_B1_translation + T_B2_translation:

# ----------------Potential energy (V)-----------------
height_of_B1 := r_B1 &o (N &ev [0, 0, 1]):
height_of_B2 := r_B2 &o (N &ev [0, 0, 1]):
V_B1 := m0*g*height_of_B1:
V_B2 := m0*g*height_of_B2:
V_S1 := (1/2)*k0*(R*q1)^2:
V_S2 := (1/2)*k1*(R*q3)^2:
V := V_B1 + V_B2 + V_S1 + V_S2:

# ------------Construct the Lagraingian----------------
L := T-V:
L := subs(q1t = u1, q2t = u2, q3t = u3, L):

# ------Construct the Euler-Lagrange equations---------
dLdq1 := diff(L, q1):
dLdu1 := diff(L, u1):
dtdLdu1 := &dt(dLdu1):
eq1 := dtdLdu1-dLdq1 = 0:
eq1 := subs(q1t = u1, q2t = u2, q3t = u3, eq1):

dLdq2 := diff(L, q2):
dLdu2 := diff(L, u2):
dtdLdu2 := &dt(dLdu2):
eq2 := dtdLdu2-dLdq2 = 0:
eq2 := subs(q1t = u1, q2t = u2, q3t = u3, eq2):

dLdq3 := diff(L, q3):
dLdu3 := diff(L, u3):
dtdLdu3 := &dt(dLdu3):
eq3 := dtdLdu3-dLdq3 = 0:
eq3 := subs(q1t = u1, q2t = u2, q3t = u3, eq3):

kde := {q1t = u1, q2t = u2, q3t = u3}:

eq1 := {u1t = solve(eq1, u1t)}:
eq2 := {u2t = solve(eq2, u2t)}:
eq3 := {u3t = solve(eq3, u3t)}:

eqs := eq1 union eq2 union eq3 union kde:

# -------------The initial conditions------------------
Initcond := {q1(0) = Pi/2-0.3,
             q2(0) = 0,
             q3(0) = 0.1,
             u1(0) = 0,
             u2(0) = 1,
             u3(0) = 0}:

eqst := subs(toTimeFunction, eqs):

# -------------------The constants---------------------
param := {R = 5,
          r = 0.6,
          H0 = Pi-0.1,
          H1 = Pi+0.1,
          l0 = 5*Pi/2,
          l1 = 5*Pi/2,
          m0 = 100,
          m1 = 10,
          k0 = 900,
          k1 = 300,
          g = 9.82}:

eqst := subs(param, eqst):

# ---------Numericall solve the system-----------------
ff := dsolve(eqst union Initcond, {q1(t),
                                   q2(t),
                                   q3(t),
                                   u1(t),
                                   u2(t),
                                   u3(t)}, type = numeric, maxfun = 0):

numberOfPoints := 1800:
integrationTime := 30:
_plot := odeplot(ff,
    [[t,q1(t)],
    [t,q2(t)],
    [t,q3(t)],
    [t,u1(t)],
    [t,u2(t)],
    [t,u3(t)]],
    0..integrationTime,numpoints=numberOfPoints):

# --------Exporting the data to a .txt file------------
# The code below simply exports the data in _plot
# in a text file name 'maple_data.txt' with the
# following format (columns are seperated by space):

# | t   | q_1 | ... | q_n | u_1 | ... | u_n |
# |-----+-----+-----+-----+-----+-----+-----|
# | t_0 |     |     |     |     |     |     |
# | t_1 |     |     |     |     |     |     |
# | .   |     |     |     |     |     |     |
# | .   |     |     |     |     |     |     |
# | t_k |     |     |     |     |     |     |

t_data := ((plottools[getdata](_plot, "points")[1])[3])[[1..numberOfPoints],1]:
q1_data := ((plottools[getdata](_plot, "points")[1])[3])[[1..numberOfPoints],2]:
q2_data := ((plottools[getdata](_plot, "points")[2])[3])[[1..numberOfPoints],2]:
q3_data := ((plottools[getdata](_plot, "points")[3])[3])[[1..numberOfPoints],2]:
u1_data := ((plottools[getdata](_plot, "points")[4])[3])[[1..numberOfPoints],2]:
u2_data := ((plottools[getdata](_plot, "points")[5])[3])[[1..numberOfPoints],2]:
u3_data := ((plottools[getdata](_plot, "points")[6])[3])[[1..numberOfPoints],2]:

_data := Concatenate(2,
                     t_data,
                     q1_data,
                     q2_data,
                     q3_data,
                     u1_data,
                     u2_data,
                     u3_data):

writedata("/path/to/maple_data.txt",convert(_data,array),float):
