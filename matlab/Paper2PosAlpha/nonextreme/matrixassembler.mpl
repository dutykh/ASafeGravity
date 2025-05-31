# Matrices assembling for Asymptotic Safe Gravity (ASG) black holes's quasi-normal modes computation
# Spin s \in {0, 1, 2} for Bonanno-Reuter black holes

MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  N::integer, # N : number of Tchebyshev modes
  s::numeric, # s : spin of perturbations
  L::numeric, # L : angular momentum, L >= s
  g::numeric, # g : gamma
  M::numeric, # M : Black-Hole mass
  p::string   # p : string of characters containing the path where we save the assembled matrices
  )
  local F::function, Fp::function, rh::numeric, alpha::numeric, epsilon::numeric, lambda::numeric, a::numeric, V::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string, T_expr::list, DT_expr::list, D2T_expr::list, T_mat::Matrix, DT_mat::Matrix, D2T_mat::Matrix, L00_vec::Vector, L01_vec::Vector, L02_vec::Vector, L10_vec::Vector, L11_vec::Vector, L12_vec::Vector, L20_vec::Vector, L21_vec::Vector, L22_vec::Vector:

  with(LinearAlgebra): # To work with matrices

  Digits := d: # Setup multiprecision arithmetics

  alpha := 118/(15*Pi):
  epsilon := 1 - s*s:
  lambda := L*(L + 1):

  # Determination of the event horizon:
  rh := fsolve(1 - 2*M*rh^2/(rh^3 + alpha*(M*g + rh)) = 0, rh = 9.03):
  printf("Found BH's event horizon rh = %f\n", rh):

  # Physical parameters:
  a := 2*M*rh/(-4*M*rh + 3*rh^2 + alpha):

  # Other functions appearing in the model:
  F := proc(y) option remember;
    1 - 8*M*rh^2*(1 - y)/(8*rh^3 + 2*M*rh^2*(1 - y)^3 + 2*alpha*rh*(1 - y)^2 + (-rh^3 - alpha*rh)*(1 - y)^3):
  end proc:

  Fp := proc(y) option remember;
    8*rh*M*((-1/2*y^3 + 3/2*y^2 - 3/2*y + 5/2)*rh^2 + M*(-1 + y)^3*rh - y*alpha*(-1 + y)^2/2)/((M*rh - 1/2*rh^2 - 1/2*alpha)*y^3 + (-3*M*rh + 3/2*rh^2 + 1/2*alpha)*y^2 + (3*M*rh - 3/2*rh^2 + 1/2*alpha)*y - M*rh - (7*rh^2)/2 - alpha/2)^2:
  end proc:

  V := proc(y) option remember;
    ((1 - y)^2)/4*F(y)*(epsilon*(1 - y)*Fp(y) + lambda):
  end proc:

  # Definition of the 2nd order ODE coefficients:
  L00 := proc(y) option remember;
    -4*V(y)/((1 - y)^2*(y + 1)):
  end proc:

  L01 := proc(y) option remember;
    (1 - y)*F(y)*((1 - y)*Fp(y) - 2*F(y))/(y + 1):
  end proc:

  L02 := proc(y) option remember;
    (1 - y)^2*F(y)^2/(y + 1):
  end proc:

  L10 := proc(y) option remember;
    F(y)*((y + 1)*(4*rh + 2*(2*M - rh)*(1 - y) - (a*rh + 2*M)*(1 - y)^2)*Fp(y) - F(y)*(8*M - (3 + y)*(1 - y)*(a*rh + 2*M)))/(y + 1)^3:
  end proc:

  L11 := proc(y) option remember;
    2*F(y)^2*(4*rh + 2*(2*M - rh)*(1 - y) - (a*rh + 2*M)*(1 - y)^2)/(y + 1)^2:
  end proc:

  L12 := proc(y) option remember;
    0:
  end proc:

  L20 := proc(y) option remember;
    4*rh^2/((y + 1)*(1 - y)^2) - F(y)^2*(2*M*(-y^2 + 1) - rh*(a*y^2 - 2*y*(1 + a) + a - 2))^2/((y + 1)^3*(1 - y)^2):
  end proc:

  L21 := proc(y) option remember;
    0:
  end proc:

  L22 := proc(y) option remember;
    0:
  end proc:

  # Precompute collocation points
  xi := [seq(evalf(cos((2.0*i-1.0)*Pi/(2.0*N))), i=1..N)]:

  # Precompute Chebyshev polynomial expressions and derivatives
  T_expr := [seq(ChebyshevT(j, x), j=0..N-1)]:
  DT_expr := [seq(diff(T_expr[j+1], x), j=0..N-1)]:
  D2T_expr := [seq(diff(T_expr[j+1], x$2), j=0..N-1)]:

  # Precompute matrices T_mat, DT_mat, D2T_mat
  T_mat := Matrix(N, N, (i,j) -> evalf(subs(x=xi[i], T_expr[j]))):
  DT_mat := Matrix(N, N, (i,j) -> evalf(subs(x=xi[i], DT_expr[j]))):
  D2T_mat := Matrix(N, N, (i,j) -> evalf(subs(x=xi[i], D2T_expr[j]))):

  # Precompute coefficient vectors
  L00_vec := Vector(N, i -> evalf(L00(xi[i]))):
  L01_vec := Vector(N, i -> evalf(L01(xi[i]))):
  L02_vec := Vector(N, i -> evalf(L02(xi[i]))):
  L10_vec := Vector(N, i -> evalf(L10(xi[i]))):
  L11_vec := Vector(N, i -> evalf(L11(xi[i]))):
  L12_vec := Vector(N, i -> evalf(L12(xi[i]))): # Will be zero
  L20_vec := Vector(N, i -> evalf(L20(xi[i]))):
  L21_vec := Vector(N, i -> evalf(L21(xi[i]))): # Will be zero
  L22_vec := Vector(N, i -> evalf(L22(xi[i]))): # Will be zero

  # Assemble matrices efficiently
  M0 := Matrix(N, N, (i,j) -> L00_vec[i]*T_mat[i,j] + L01_vec[i]*DT_mat[i,j] + L02_vec[i]*D2T_mat[i,j]):
  M1 := Matrix(N, N, (i,j) -> L10_vec[i]*T_mat[i,j] + L11_vec[i]*DT_mat[i,j] + L12_vec[i]*D2T_mat[i,j]):
  M2 := Matrix(N, N, (i,j) -> L20_vec[i]*T_mat[i,j] + L21_vec[i]*DT_mat[i,j] + L22_vec[i]*D2T_mat[i,j]):

  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(N, string):

  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):

end proc:
