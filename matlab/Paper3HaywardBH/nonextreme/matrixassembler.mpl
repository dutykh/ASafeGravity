# Matrices assembling for Asymptotic Safe Gravity (ASG) black holes's quasi-normal modes computation
# Spin s \in {0, 1, 2} for Hayward BHs (non-extreme case)

MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  N::integer, # N : number of Tchebyshev modes
  s::numeric, # s : spin of perturbations
  L::numeric, # L : angular momentum, L >= s
  g::numeric, # g : gamma
  p::string   # p : string of characters containing the path where we save the assembled matrices
  )
  local F::function, Fp::function, rhoh::numeric, alpha::numeric, epsilon::numeric, lambda::numeric, a:: numeric, V::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi, path::string, nstr::string, get_rhoh_guess::procedure, T_expr, DT_expr, D2T_expr, T_mat, DT_mat, D2T_mat, L00_vec, L01_vec, L02_vec, L10_vec, L11_vec, L12_vec, L20_vec, L21_vec, L22_vec, all_roots_list, real_roots, guess, closest_root, min_dist, dist:

  with(LinearAlgebra): # To work with matrices

  Digits := d: # Setup multiprecision arithmetics

  epsilon := 1 - s*s:
  lambda := L*(L + 1):

  # Procedure to get an initial guess for rhoh from a lookup table
  get_rhoh_guess := proc(g_val::numeric)
    local g_rhoh_table, g_vals, i;
    g_rhoh_table := [
      [0, 1.99], [0.02, 1.995], [0.04, 1.990], [0.06, 1.985], [0.08, 1.980], [0.10, 1.974],
      [4/27, 1.961], [8/27, 1.920], [0.30, 1.918], [12/27, 1.873], [0.50, 1.855],
      [16/27, 1.821], [0.70, 1.779], [20/27, 1.761], [24/27, 1.688], [0.90, 1.682], [1.0, 1.618],
      [28/27, 1.590], [1.10, 1.530], [1.18, 1.384], [32/27, 1.33]
    ];
    g_vals := map(x -> x[1], g_rhoh_table);
    for i from 1 to numelems(g_vals) do
      if g_val <= g_vals[i] then
        return g_rhoh_table[i][2];
      end if;
    end do;
    return g_rhoh_table[numelems(g_rhoh_table)][2]; # Fallback to the last value
  end proc:

  # Determination of the event horizon:
  all_roots_list := [allvalues(RootOf(_Z^3 - 2*_Z^2 + g, _Z))]:
  real_roots := select(type, all_roots_list, realcons):

  if numelems(real_roots) = 0 then
    error "No real roots found for the event horizon equation for g=%a", g;
  elif numelems(real_roots) = 1 then
    rhoh := real_roots[1]:
  else
    guess := get_rhoh_guess(g):
    # Find the root closest to the guess
    closest_root := real_roots[1]:
    min_dist := abs(closest_root - guess):
    for i from 2 to numelems(real_roots) do
        dist := abs(real_roots[i] - guess):
        if dist < min_dist then
            min_dist := dist:
            closest_root := real_roots[i]:
        end if:
    end do:
    rhoh := closest_root:
  end if:
  printf("Found BH's event horizon rh = %f\n", rhoh):

  # Physical parameters:
  a := 2/(3*rhoh - 4):

  # Other functions appearing in the model:
  F := proc(y) option remember;
    (1 + y)*((rhoh - 2)*y^2 + (-4*rhoh + 8)*y + 7*rhoh - 6)/((rhoh - 2)*y^3 + (-3*rhoh + 6)*y^2 + (3*rhoh - 6)*y + 7*rhoh + 2):
  end proc:

  Fp := proc(y) option remember;
    ((-16*rhoh + 32)*y^3 + (48*rhoh - 96)*y^2 + (-48*rhoh + 96)*y + 80*rhoh - 32)/((rhoh - 2)*y^3 + (-3*rhoh + 6)*y^2 + (3*rhoh - 6)*y + 7*rhoh + 2)^2:
  end proc:

  V := proc(y) option remember;
    ((1 - y)^2)/4*F(y)*(epsilon*(1 - y)*Fp(y) + lambda):
  end proc:

  # Definition of the 2nd order ODE coefficients:
  L00 := proc(y) option remember;
    -4*V(y)/((1 - y)^2*(1 + y)):
  end proc:

  L01 := proc(y) option remember;
    (1 - y)*F(y)*((1 - y)*Fp(y) - 2*F(y))/(1 + y):
  end proc:

  L02 := proc(y) option remember;
    (1 - y)^2*(F(y)^2)/(1 + y):
  end proc:

  L10 := proc(y) option remember;
    F(y)*((1 + y)*(4*rhoh + 2*(2 - rhoh)*(1 - y) - (a*rhoh + 2)*(1 - y)^2)*Fp(y) - F(y)*(8 - (3 + y)*(1 - y)*(a*rhoh + 2)))/(1 + y)^3:
  end proc:

  L11 := proc(y) option remember;
    2*F(y)^2*(4*rhoh + 2*(2 - rhoh)*(1 - y) - (a*rhoh + 2)*(1 - y)^2)/(1 + y)^2:
  end proc:

  L12 := proc(y) option remember;
    0:
  end proc:

  L20 := proc(y) option remember;
    4*rhoh^2/((1 + y)*(1 - y)^2) - F(y)^2*(2*(-y^2 + 1) - rhoh*(a*y^2 - 2*y*(1 + a) + a - 2))^2/((1 + y)^3*(1 - y)^2):
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