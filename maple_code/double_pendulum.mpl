# By Filip Strand and Jakob Arnoldsson
# May 2016
# (used with Maple 2015 and Sophia for Maple R6 - 13 june 2000) 
#
read "/path/to/SophiaV6.txt":
read "/path/to/Graphics.txt":
with(plots): with(plottools): with(linalg): with(ArrayTools):

# --------------Declaring time depencece---------------
dependsTime(q1, q2, q3, q4, u1, u2, u3, u4):

# --------------Coordinate Systems---------------------
&rot([N, A, 2, q1]):
&rot([A, B, 2, q2]):

# ------------------Position vectors-------------------
r_O := (N &ev [0,0,0]):
r_A := (N &ev [l0+q4,0,0]):
r_B := (N &ev [l0+q4,0,0]) &++ (A &ev [0, 0, -l1/2]):
r_C := (N &ev [l0+q4,0,0]) &++ (A &ev [0, 0, -l1]):
r_D := (N &ev [l0+q4,0,0]) &++ (A &ev [0, 0, -l1]) &++ (B &ev [0, 0, -(l2 + q3)]):

# ----------------Velocity vectors --------------------
v_B := simplify(N &fdt r_B):
v_D := simplify(N &fdt r_D):

# ---------------Angular momentum (H)------------------
omega := <0,q1t,0>:
I_rod := << (1/12)*m1*l1^2 | 0 | 0 >,
          < 0 | (1/12)*m1*l1^2 | 0 >,
          < 0 | 0 | 0 >>:
H_rod := I_rod.omega:

# ----------------Kinetic energy (T)-------------------
T_rod_rotation := (1/2)*(omega^+).H_rod:
T_rod_translation := (1/2)*m1*(v_B &o v_B):
T_ball_translation := (1/2)*m2*(v_D &o v_D):
T := T_rod_rotation + T_rod_translation + T_ball_translation:

# ----------------Potential energy (V)-----------------
height_of_r_B := r_B &o (N &ev [0, 0, 1]):
height_of_r_D := r_D &o (N &ev [0, 0, 1]):
V_rod := m1*g*height_of_r_B:
V_ball := m2*g*height_of_r_D:
V_S1 := (1/2)*k1*q4^2:
V_S2 := (1/2)*k2*q3^2:
V := V_rod + V_ball + V_S1 + V_S2:

# ------------Construct the Lagraingian----------------
L := T-V:
L := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, L):

# ------Construct the Euler-Lagrange equations---------
dLdq1 := diff(L, q1):
dLdu1 := diff(L, u1):
dtdLdu1 := &dt(dLdu1):
eq1 := dtdLdu1-dLdq1 = 0:
eq1 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq1):

dLdq2 := diff(L, q2):
dLdu2 := diff(L, u2):
dtdLdu2 := &dt(dLdu2):
eq2 := dtdLdu2-dLdq2 = 0:
eq2 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq2):

dLdq3 := diff(L, q3):
dLdu3 := diff(L, u3):
dtdLdu3 := &dt(dLdu3):
eq3 := dtdLdu3-dLdq3 = 0:
eq3 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq3):

dLdq4 := diff(L, q4):
dLdu4 := diff(L, u4):
dtdLdu4 := &dt(dLdu4):
eq4 := dtdLdu4-dLdq4 = 0:
eq4 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq4):

kde := {q1t = u1, q2t = u2, q3t = u3, q4t = u4}:

eq1 := {u1t = solve(eq1, u1t)}:
eq2 := {u2t = solve(eq2, u2t)}:
eq3 := {u3t = solve(eq3, u3t)}:
eq4 := {u4t = solve(eq4, u4t)}:

eqs := eq1 union eq2 union eq3 union eq4 union kde:

# -------------The initial conditions------------------
Initcond := {q1(0) = -Pi/4,
             q2(0) = 0,
             q3(0) = 0,
             q4(0) = 0,
             u1(0) = 0,
             u2(0) = 0,
             u3(0) = 0,
             u4(0) = 0}:

eqst := subs(toTimeFunction, eqs):

# -------------------The constants---------------------
param := {l0 = 3,
          l1 = 7,
          l2 = 3,
          m1 = 1,
          m2 = 0.5,
          k1 = 50,
          k2 = 50,
          g = 9.82}:

eqst := subs(param, eqst):

# ---------Numericall solve the system-----------------
ff := dsolve(eqst union Initcond, {q1(t),
                                   q2(t),
                                   q3(t),
                                   q4(t),
                                   u1(t),
                                   u2(t),
                                   u3(t),
                                   u4(t)}, type = numeric, maxfun = 0):

numberOfPoints:= 720:
integrationTime:= 30:
_plot := odeplot(ff,
                 [[t,q1(t)],
                  [t,q2(t)],
                  [t,q3(t)],
                  [t,q4(t)],
                  [t,u1(t)],
                  [t,u2(t)],
                  [t,u3(t)],
                  [t,u4(t)]],0..integrationTime,numpoints=numberOfPoints):

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
q4_data := ((plottools[getdata](_plot, "points")[4])[3])[[1..numberOfPoints],2]:
u1_data := ((plottools[getdata](_plot, "points")[5])[3])[[1..numberOfPoints],2]:
u2_data := ((plottools[getdata](_plot, "points")[6])[3])[[1..numberOfPoints],2]:
u3_data := ((plottools[getdata](_plot, "points")[7])[3])[[1..numberOfPoints],2]:
u4_data := ((plottools[getdata](_plot, "points")[8])[3])[[1..numberOfPoints],2]:

_data := Concatenate(2,
                     t_data,
                     q1_data,
                     q2_data,
                     q3_data,
                     q4_data,
                     u1_data,
                     u2_data,
                     u3_data,
                     u4_data):

writedata("/path/to/maple_data.txt",convert(_data,array),float):
