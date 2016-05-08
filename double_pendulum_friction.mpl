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

# ---------Construct the Generalized forces------------
F1 := N &ev [-k1*q4,0,0]:
F2 := N &ev [0,0,-m1*g]:
F3 := N &to (B &ev [0,0,-k2*q3]):
F4 := N &to (B &ev [0,0,k2*q3]):
F5 := N &ev [0,0,-m2*g]:

Q1 := F1 &o [map(diff,r_A[1],q1),r_A[2]] +
      F2 &o [map(diff,r_B[1],q1),r_B[2]] +
      F3 &o [map(diff,r_C[1],q1),r_C[2]] +
      F4 &o [map(diff,r_D[1],q1),r_D[2]] +
      F5 &o [map(diff,r_D[1],q1),r_D[2]] -
      C1*q1t:

Q2 := F1 &o [map(diff,r_A[1],q2),r_A[2]] +
      F2 &o [map(diff,r_B[1],q2),r_B[2]] +
      F3 &o [map(diff,r_C[1],q2),r_C[2]] +
      F4 &o [map(diff,r_D[1],q2),r_D[2]] +
      F5 &o [map(diff,r_D[1],q2),r_D[2]] -
      C2*q2t:

Q3 := F1 &o [map(diff,r_A[1],q3),r_A[2]] +
      F2 &o [map(diff,r_B[1],q3),r_B[2]] +
      F3 &o [map(diff,r_C[1],q3),r_C[2]] +
      F4 &o [map(diff,r_D[1],q3),r_D[2]] +
      F5 &o [map(diff,r_D[1],q3),r_D[2]]:

Q4 := F1 &o [map(diff,r_A[1],q4),r_A[2]] +
      F2 &o [map(diff,r_B[1],q4),r_B[2]] +
      F3 &o [map(diff,r_C[1],q4),r_C[2]] +
      F4 &o [map(diff,r_D[1],q4),r_D[2]] +
      F5 &o [map(diff,r_D[1],q4),r_D[2]]:

# ------Construct the Euler-Lagrange equations---------
T := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, T):

dTdq1 := diff(T, q1):
dTdu1 := diff(T, u1):
dtdTdu1 := &dt(dTdu1):
eq1 := dtdTdu1-dTdq1-Q1 = 0:
eq1 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq1):

dTdq2 := diff(T, q2):
dTdu2 := diff(T, u2):
dtdTdu2 := &dt(dTdu2):
eq2 := dtdTdu2-dTdq2-Q2 = 0:
eq2 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq2):

dTdq3 := diff(T, q3):
dTdu3 := diff(T, u3):
dtdTdu3 := &dt(dTdu3):
eq3 := dtdTdu3-dTdq3-Q3 = 0:
eq3 := subs(q1t = u1, q2t = u2, q3t = u3, q4t = u4, eq3):

dTdq4 := diff(T, q4):
dTdu4 := diff(T, u4):
dtdTdu4 := &dt(dTdu4):
eq4 := dtdTdu4-dTdq4-Q4 = 0:
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
          C1 = 10,
          C2 = 10,
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

numberOfPoints:= 1800:
integrationTime:= 60:
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
