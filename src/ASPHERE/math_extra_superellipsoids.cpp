// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Jacopo Bilotto (EPFL), Jibril Coulibaly (??)
------------------------------------------------------------------------- */

#include "math_extra_superellipsoids.h"
#include "math_extra.h"
#include <cmath>

// #include "math_special.h"
// #include "math_const.h"

// #include <algorithm>
// #include <cstdio>
// #include <cstring>

extern "C" { // General Matrices
    void dgetrf_(const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info); // Factorize
    void dgetrs_(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info); // Solve (using factorzation)
}

namespace MathExtraSuperellipsoids {

static constexpr int ITERMAX_NEWTON = 100;
static constexpr double CONVERGENCE_NEWTON = 1e-6;
static constexpr int ITERMAX_LINESEARCH = 10;
static constexpr double PARAMETER_LINESEARCH = 1e-4;
static constexpr double CUTBACK_LINESEARCH = 0.5;
static constexpr double CONVERGENCE_OVERLAP = 1e-8;
static constexpr unsigned int ITERMAX_OVERLAP = 20;
static constexpr double MINSLOPE_OVERLAP = 1e-12;

/* ----------------------------------------------------------------------
   beta function B(x,y) = Gamma(x) * Gamma(y) / Gamma(x+y)
------------------------------------------------------------------------- */
double beta_func(double a, double b) {
    return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}
// TODO: the reason why I had codded the beta function from scratch is because LAMMPS must be guaranteed to work with some older standard of C++
//       I don't remember which one exactly (C++14 I think) but this standard does not have gamma() of beta() in the <cmath> implementation
//       TBD if the code above will be accepted or if we need to fall back to the implementation of beta I copied from Cephes
/* ----------------------------------------------------------------------
   Volume of superellipsoid
   source https://cse.buffalo.edu/~jryde/cse673/files/superquadrics.pdf
------------------------------------------------------------------------- */

void volume_superellipsoid(const double *blockiness, const double *shape, double volume)
{
  const double eps1 = 2.0 / blockiness[0]; // shape exponent in latitude direction
  const double eps2 = 2.0 / blockiness[1]; // shape exponent in longitude direction
  volume = 2.0*shape[0]*shape[1]*shape[2]*eps1*eps2*
      beta_func(0.5*eps1, eps1 + 1.0)*
      beta_func(0.5*eps2, 0.5*eps2 + 1.0);
}

/* ----------------------------------------------------------------------
   inertia tensor of superellipsoid
   source https://cse.buffalo.edu/~jryde/cse673/files/superquadrics.pdf
------------------------------------------------------------------------- */
void inertia_superellipsoid(const double *shape, const double *blockiness, double density, double *inertia)

{

  const double eps1 = 2.0 / blockiness[0]; // shape exponent in latitude direction
  const double eps2 = 2.0 / blockiness[1]; // shape exponent in longitude direction

  const double a1 = shape[0];
  const double a2 = shape[1];
  const double a3 = shape[2];
  const double I_xx = 0.5*a1*a2*a3*eps1*eps2*(a2*a2*beta_func(1.5*eps2, 0.5*eps2)*beta_func(0.5*eps1, 2.0*eps1+1.0)+
      4.0*a3*a3*beta_func(0.5*eps2, 0.5*eps2+1.0)*beta_func(1.5*eps1, eps1+1.0)) * density;
  const double I_yy = 0.5*a1*a2*a3*eps1*eps2*(a1*a1*beta_func(1.5*eps2, 0.5*eps2)*beta_func(0.5*eps1, 2.0*eps1+1.0)+
      4.0*a3*a3*beta_func(0.5*eps2, 0.5*eps2+1.0)*beta_func(1.5*eps1, eps1+1.0)) * density;
  const double I_zz = 0.5*a1*a2*a3*eps1*eps2*(a1*a1 + a2*a2)*
      beta_func(1.5*eps2, 0.5*eps2)*beta_func(0.5*eps1, 2.0*eps1+1.0) * density;

  inertia[0] = I_xx;
  inertia[1] = I_yy;
  inertia[2] = I_zz;
}


/* ----------------------------------------------------------------------
   curvature of superellipsoid
   source https://en.wikipedia.org/wiki/Mean_curvature
------------------------------------------------------------------------- */

void mean_curvature_superellipsoid(const double *shape, const double *blockiness, const double* quat, const double *global_point, double curvature)
{
  // this code computes the mean curvature on the superellipsoid surface
  // for the given global point
  double local_point[3],hessian[3][3], nablaF[3], f, normal[3];
  global2local_vector(global_point, quat, local_point); 
  shape_function_local(shape, blockiness, quat, local_point, f);
  double koef = pow(fabs(0.5), std::max(blockiness[0], blockiness[1])-2.0);
  double alpha = 1.0 / pow(fabs(f/koef + 1.0), 1.0/blockiness[0]);
  for(int i = 0; i < 3; i++)
    local_point[i] *= alpha;
  shape_function_local_grad(shape, blockiness, quat, local_point, nablaF);
  shape_function_local_hessian(shape, blockiness, quat, local_point, hessian);
  MathExtra::normalize3(nablaF, normal);
  double temp[3];
  MathExtra::matvec(hessian, normal, temp);
  double F_mag = sqrt(MathExtra::dot3(nablaF, nablaF));
  curvature = fabs(MathExtra::dot3(normal, temp) - (hessian[0][0] + hessian[1][1] + hessian[2][2])) / fabs(2.0 * F_mag);
}

void gaussian_curvature_superellipsoid(const double *shape, const double *blockiness, const double* quat, const double *global_point, double curvature)
{
  // this code computes the gaussian curvature coefficient
  // for the given global point
  double local_point[3],hessian[3][3], nablaF[3], f, normal[3];
  global2local_vector(global_point, quat, local_point); 
  shape_function_local(shape, blockiness, quat, local_point, f);
  double koef = pow(fabs(0.5), std::max(blockiness[0], blockiness[1])-2.0);
  double alpha = 1.0 / pow(fabs(f/koef + 1.0), 1.0/blockiness[0]);
  for(int i = 0; i < 3; i++)
    local_point[i] *= alpha;
  shape_function_local_grad(shape, blockiness, quat, local_point, nablaF);
  shape_function_local_hessian(shape, blockiness, quat, local_point, hessian);
  MathExtra::normalize3(nablaF, normal);
  double temp[3];
  MathExtra::matvec(hessian, normal, temp);
  double F_mag = sqrt(MathExtra::dot3(nablaF, nablaF));

  double fx = nablaF[0];
  double fy = nablaF[1];
  double fz = nablaF[2];

  double fxx = hessian[0][0];
  double fxy = hessian[0][1];
  double fxz = hessian[0][2];

  double fyy = hessian[1][1];
  double fyz = hessian[1][2];

  double fzz = hessian[2][2];

  double mat[4][4] = {
    {fxx, fxy, fxz, fx},
    {fxy, fyy, fyz, fy},
    {fxz, fyz, fzz, fz},
    {fx,  fy,  fz, 0.0} 
  };

    double K = -det4_M44_zero(mat) / (F_mag*F_mag*F_mag*F_mag);
    curvature =  sqrt(fabs(K));
}


/* ----------------------------------------------------------------------
   express local (particle level) to global (system level) coordinates
------------------------------------------------------------------------- */

void local2global_vector(const double v[3], const double *quat, double global_v[3]){

   MathExtra::quatrotvec(const_cast<double*>(quat) , const_cast<double*>(v), global_v);
};

void local2global_matrix(const double m[3][3], const double *quat, double global_m[3][3]){
    double rot[3][3],  temp[3][3];
    MathExtra::quat_to_mat(const_cast<double*>(quat), rot);
    MathExtra::times3(rot, m, temp);
    MathExtra::transpose_times3(rot, temp, global_m);
};

  
/* ----------------------------------------------------------------------
   express global (system level) to local (particle level) coordinates
------------------------------------------------------------------------- */

void global2local_vector(const double *v, const double *quat, double *local_v){

    double qc[4];
    MathExtra::qconjugate(const_cast<double*>(quat), qc);
    MathExtra::quatrotvec(qc, const_cast<double*>(v), local_v);

};


void global2local_matrix(const double m[3][3], const double *quat, double local_m[3][3]){
    double rot[3][3], temp[3][3];
    MathExtra::quat_to_mat(quat, rot);
    MathExtra::transpose_times3(rot, m, temp);
    MathExtra::times3(temp, rot, local_m);
}

/* ----------------------------------------------------------------------
   shape function computations for superellipsoids
------------------------------------------------------------------------- */

void shape_function_local(const double *shape, const double *block, const double *quat, const double *point, double local_f){
  const double n1 = block[0], n2 = block[1];
  
  local_f = pow( pow(abs(point[0]/shape[0]), n2) + pow(abs(point[1]/shape[1]), n2) , n1/ n2) + pow(abs(point[2]/shape[2]), n1)  - 1.0;
};

void shape_function_global(const double *shape, const double *block, const double *quat, const double *point, double global_f){
  double local_point[3];
  global2local_vector(const_cast<double*>(point), const_cast<double*>(quat), local_point);
  shape_function_local(shape, block, quat, local_point, global_f);
};

void shape_function_local_grad(const double *shape, const double *block, const double *quat, const double *point, double *local_grad){
  // point is in local coordinates
  const double n1 = block[0], n2 = block[1];
  const double ainv = 1.0 / shape[0];
  const double binv = 1.0 / shape[1];
  const double cinv = 1.0 / shape[2];

  const double nu = pow(abs(point[0] * ainv), n2) + pow(abs(point[1] * binv), n2);
  const double nu_12 = pow(nu, n1 / n2 - 1.0);

  local_grad[0] = n1*ainv * pow(abs(point[0] * ainv), n2 - 1.0) * nu_12 * copysign(1.0, point[0]);
  local_grad[1] = n1*binv * pow(abs(point[1] * binv), n2 - 1.0) * nu_12 * copysign(1.0, point[1]);
  local_grad[2] = n1*cinv * pow(abs(point[2] * cinv), n1 - 1.0) * copysign(1.0, point[2]);

};

void shape_function_local_hessian(
  const double *shape, const double *block, const double *quat, const double *point, double local_hess[3][3]) {
  const double n1 = block[0], n2 = block[1];
  const double ainv = 1.0 / shape[0];
  const double binv = 1.0 / shape[1];
  const double cinv = 1.0 / shape[2];

  const double nu = pow(abs(point[0] * ainv), n2) + pow(abs(point[1] * binv), n2);
  const double nu_12_1 = pow(nu, n1 / n2 - 1.0);
  const double nu_12_2 = pow(nu, n1 / n2 - 2.0);

  local_hess[0][2] = local_hess[2][0] = local_hess[1][2] = local_hess[2][1] =0;

  local_hess[0][0] = n1 * (n2 - 1) * ainv * ainv * pow(abs(point[0] * ainv), n2 - 2.0)* nu_12_1 +
                     n1 * (n1 - n2) * ainv * ainv * pow(abs(point[0] * ainv), 2*n2 - 2.0)* nu_12_2;

  local_hess[1][1] = n1 * (n2 - 1) * binv * binv * pow(abs(point[1] * binv), n2 - 2.0)* nu_12_1 +
                     n1 * (n1 - n2) * ainv * ainv * pow(abs(point[1] * binv), 2*n2 - 2.0)* nu_12_2;

  local_hess[2][2] = n1 * (n1 - 1) * cinv * cinv * pow(abs(point[2] * cinv), n1-2);

  local_hess[0][1] = n1 * (n1 - n2) * ainv * binv * pow(abs(point[0]*ainv), n2 - 1) *
                     pow(abs(point[1]*binv), n2 -1) * pow(nu, n1 / n2 - 2) * copysign(1.0, shape[0] * shape[1]); 
                
  }


double shape_and_derivatives_local(const double* xlocal, const double* shape, const double* block, const int flag, double* grad, double hess[3][3]) {
  double shapefunc;
  switch (flag) {
    case 0: {
      shapefunc = shape_and_derivatives_local_ellipsoid(xlocal, shape, grad, hess);
      break;
    }
    case 1: {
      shapefunc = shape_and_derivatives_local_n1equaln2(xlocal, shape, block[0], grad, hess);
      break;
    }
    case 2: {
      shapefunc = shape_and_derivatives_local_superquad(xlocal, shape, block, grad, hess);
      break;
    }
  }
  return shapefunc;
}

// General case for n1 != n2 > 2
double shape_and_derivatives_local_superquad(const double* xlocal, const double* shape, const double* block, double* grad, double hess[3][3]) {
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = std::fabs(xlocal[0] * a_inv);
  double y_b = std::fabs(xlocal[1] * b_inv);
  double z_c = std::fabs(xlocal[2] * c_inv);
  double n1 = block[0];
  double n2 = block[1];
  double x_a_pow_n2_m2 = std::pow(x_a, n2 - 2.0);
  double x_a_pow_n2_m1 = x_a_pow_n2_m2 * x_a;
  double y_b_pow_n2_m2 = std::pow(y_b, n2 - 2.0);
  double y_b_pow_n2_m1 = y_b_pow_n2_m2 * y_b;

  double nu = (x_a_pow_n2_m1 * x_a) + (y_b_pow_n2_m1 * y_b);
  double nu_pow_n1_n2_m2 = std::pow(nu, n1/n2 - 2.0);
  double nu_pow_n1_n2_m1 = nu_pow_n1_n2_m2 * nu;

  double z_c_pow_n1_m2 = std::pow(z_c, n1 -2.0);
  double z_c_pow_n1_m1 = z_c_pow_n1_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n1 * a_inv * x_a_pow_n2_m1 * nu_pow_n1_n2_m1 * signx;
  grad[1] = n1 * b_inv * y_b_pow_n2_m1 * nu_pow_n1_n2_m1 * signy;
  grad[2] = n1 * c_inv * z_c_pow_n1_m1 * signz;

  // Equation (15)
  double signxy = signx * signy;
  hess[0][0] = a_inv * a_inv * (n1 * (n2 - 1.0) * x_a_pow_n2_m2 * nu_pow_n1_n2_m1 +
                                (n1 - n2) * n1 * (x_a_pow_n2_m1 * x_a_pow_n2_m1) * nu_pow_n1_n2_m2);
  hess[1][1] = b_inv * b_inv * (n1 * (n2 - 1.0) * y_b_pow_n2_m2 * nu_pow_n1_n2_m1 +
                                (n1 - n2) * n1 * (y_b_pow_n2_m1 * y_b_pow_n2_m1) * nu_pow_n1_n2_m2);
  hess[0][1] = hess[1][0] = a_inv * b_inv * (n1 - n2) * n1 * x_a_pow_n2_m1 * y_b_pow_n2_m1 * nu_pow_n1_n2_m2 * signxy;
  hess[2][2] = c_inv * c_inv * n1 * (n1 - 1.0) * z_c_pow_n1_m2;
  hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;

  return (nu_pow_n1_n2_m1 * nu) + (z_c_pow_n1_m1 * z_c) - 1.0;
}

// Special case for n2 = n2 = n > 2
double shape_and_derivatives_local_n1equaln2(const double* xlocal, const double* shape, const double n, double* grad, double hess[3][3]) {
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = std::fabs(xlocal[0] * a_inv);
  double y_b = std::fabs(xlocal[1] * b_inv);
  double z_c = std::fabs(xlocal[2] * c_inv);
  double x_a_pow_n_m2 = std::pow(x_a, n - 2.0);
  double x_a_pow_n_m1 = x_a_pow_n_m2 * x_a;
  double y_b_pow_n_m2 = std::pow(y_b, n - 2.0);
  double y_b_pow_n_m1 = y_b_pow_n_m2 * y_b;
  double z_c_pow_n_m2 = std::pow(z_c, n - 2.0);
  double z_c_pow_n_m1 = z_c_pow_n_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n * a_inv * x_a_pow_n_m1 * signx;
  grad[1] = n * b_inv * y_b_pow_n_m1 * signy;
  grad[2] = n * c_inv * z_c_pow_n_m1 * signz;

  // Equation (15)
  double signxy = signx * signy;
  hess[0][0] = a_inv * a_inv * n * (n - 1.0) * x_a_pow_n_m2;
  hess[1][1] = b_inv * b_inv * n * (n - 1.0) * y_b_pow_n_m2;
  hess[2][2] = c_inv * c_inv * n * (n - 1.0) * z_c_pow_n_m2;
  hess[0][1] = hess[1][0] = hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;

  return (x_a_pow_n_m1 * x_a) + (y_b_pow_n_m1 * y_b) + (z_c_pow_n_m1 * z_c) - 1.0;
}


// Special case for n1 = n2 = 2
double shape_and_derivatives_local_ellipsoid(const double* xlocal, const double* shape, double* grad, double hess[3][3]) {
  double a = 2.0 / (shape[0] * shape[0]);
  double b = 2.0 / (shape[1] * shape[1]);
  double c = 2.0 / (shape[2] * shape[2]);

  // Equation (14) simplified for n1 = n2 = 2
  grad[0] = a * xlocal[0];
  grad[1] = b * xlocal[1];
  grad[2] = c * xlocal[2];

  // Equation (15)
  hess[0][0] = a;
  hess[1][1] = b;
  hess[2][2] = c;
  hess[0][1] = hess[1][0] = hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;

  return 0.5 * (grad[0]*xlocal[0] + grad[1]*xlocal[1] + grad[2]*xlocal[2]) - 1.0;
}

double shape_and_derivatives_global(const double* xc, const double R[3][3], const double* shape, const double* block, const int flag, const double* X0, double* grad, double hess[3][3]) {
  double shapefunc, xlocal[3], tmp_v[3], tmp_m[3][3];
  MathExtra::sub3(X0, xc, tmp_v);
  MathExtra::transpose_matvec(R, tmp_v, xlocal);
  shape_and_derivatives_local(xlocal, shape, block, flag, tmp_v, hess);
  MathExtra::matvec(R, tmp_v, grad);
  MathExtra::times3_transpose(hess, R, tmp_m);
  MathExtra::times3(R, tmp_m, hess);
  return shapefunc;
}

double compute_residual(const double shapefunci, const double* gradi_global, const double shapefuncj, const double* gradj_global, const double mu2, double* residual) {
  // Equation (23)
  MathExtra::scaleadd3(mu2, gradj_global, gradi_global, residual);
  residual[3] = shapefunci - shapefuncj;
  return residual[0]*residual[0] + residual[1]*residual[1] + residual[2]*residual[2] + residual[3]*residual[3];
}

void compute_jacobian(const double* gradi_global, const double hessi_global[3][3], const double* gradj_global, const double hessj_global[3][3], const double mu2, double* jacobian) {
  // Jacobian (derivative of residual)
  // 1D column-major matrix for LAPACK/linalg compatibility
  for (int row = 0 ; row < 3 ; row++) {
    for (int col = 0 ; col < 3 ; col++) {
      jacobian[row + col*4] = hessi_global[row][col] + mu2 * hessj_global[row][col];
    }
    jacobian[row + 3*4] = gradj_global[row];
  }
  for (int col = 0 ; col < 3 ; col++) {
    jacobian[3 + col*4] = gradi_global[col] - gradj_global[col];
  }
  jacobian[15] = 0.0;
}

double compute_residual_and_jacobian(const double* xci, const double Ri[3][3], const double* shapei, const double* blocki, const int flagi,
                                     const double* xcj, const double Rj[3][3], const double* shapej, const double* blockj, const int flagj,
                                     const double* X, double* shapefunc, double* residual, double* jacobian) {
  double gradi[3], hessi[3][3], gradj[3], hessj[3][3];
  shapefunc[0] = shape_and_derivatives_global(xci, Ri, shapei, blocki, flagi, X, gradi, hessi);
  shapefunc[1] = shape_and_derivatives_global(xcj, Rj, shapej, blockj, flagj, X, gradj, hessj);
  compute_jacobian(gradi, hessi, gradj, hessj, X[3], jacobian);
  return compute_residual(shapefunc[0], gradi, shapefunc[1], gradj, X[3], residual);
}


int determine_contact_point(const double* xci, const double Ri[3][3], const double* shapei, const double* blocki, const int flagi,
                            const double* xcj, const double Rj[3][3], const double* shapej, const double* blockj, const int flagj,
                            double* X0, double* nij) {
  double norm, norm_ini, shapefunc[2], residual[4], jacobian[16];
  bool converged(false);

  norm = compute_residual_and_jacobian(xci, Ri, shapei, blocki, flagi, xcj, Rj, shapej, blockj, flagj, X0, shapefunc, residual, jacobian);
  for (int iter = 0 ; iter < ITERMAX_NEWTON ; iter++) {
    norm_ini = norm;

    // Solve Newton step
    int lapack_error, ipiv[16];
    const int n = 4;
    const char trans = 'N';
    const int nrhs = 1;
    double rhs[4] = {-residual[0], -residual[1], -residual[2], -residual[3]};
    dgetrf_(&n, &n, jacobian, &n, ipiv, &lapack_error);
    if (lapack_error)
      return lapack_error;
    dgetrs_(&trans, &n, &nrhs, jacobian, &n, ipiv, rhs, &n, &lapack_error);
    if (lapack_error)
      return lapack_error;

    // Backtracking line search
    double a(1.0), X_line[4];
    int iter_ls;
    for (iter_ls = 0 ; iter_ls < ITERMAX_LINESEARCH ; iter_ls++) {
      X_line[0] = X0[0] + a * rhs[0];
      X_line[1] = X0[1] + a * rhs[1];
      X_line[2] = X0[2] + a * rhs[2];
      X_line[3] = X0[3] + a * rhs[3];

      // Line search iterates not selected for the next Newton iteration
      // do not need to compute the expensive Jacobian, only the residual.
      // We want to avoid calling `compute_residual_and_jacobian()` for each
      // line search iterate.
      // However, many intermediate variables that are costly to compute
      // are shared by the local gradient and local hessian calculations.
      // We want to avoid calling `compute_residual()` followed by `compute_jacobian()`
      // for the iterates that satisfy the descent condition.
      // To do so, we duplicate `compute_residual_and_jacobian()`, but only
      // build the global hessians if the descent condition is satisfied and
      // the iterate will be used in the next Newton step.
      // This leads to some code duplication, and still computes
      // the local hessians even when they are not necessary.
      // This seems to be an acceptable in-between of performance and clean code.
      // As most of the cost in the Hessian is in the 2 matrix products to
      // Compute the global matrix from the local one

      // One alternative would be to store the intermediate variables from
      // the local gradient calculation when calling `shape_and_gradient_local()`,
      // and re-use them during the local hessian calculation (function that 
      // calculates only the Hessian from these intermediate values would need
      // to be implemented).
      // This seems a bit clunky just to save the few multiplications of the
      // local hessian calculation, that is why I did not do it. I am open to
      // other ideas and solutions.
      // Even then, we would have some code duplication with `compute_residual_and_jacobian()`
      // So maybe I am overthinking this...

      double xilocal[3], gradi[3], hessi[3][3], xjlocal[3], gradj[3], hessj[3][3], tmp_v[3];

      MathExtra::sub3(X_line, xci, tmp_v);
      MathExtra::transpose_matvec(Ri, tmp_v, xilocal);
      shapefunc[0] = shape_and_derivatives_local(xilocal, shapei, blocki, flagi, tmp_v, hessi);
      MathExtra::matvec(Ri, tmp_v, gradi);

      MathExtra::sub3(X_line, xcj, tmp_v);
      MathExtra::transpose_matvec(Rj, tmp_v, xjlocal);
      shapefunc[1] = shape_and_derivatives_local(xjlocal, shapej, blockj, flagj, tmp_v, hessj);
      MathExtra::matvec(Rj, tmp_v, gradj);

      norm = compute_residual(shapefunc[0], gradi, shapefunc[1], gradj, X_line[3], residual);

      if (norm > norm_ini - PARAMETER_LINESEARCH * a * norm_ini) { // Armijo - Goldstein condition not met
        a *= CUTBACK_LINESEARCH;
      } else {
        X0[0] = X_line[0];
        X0[1] = X_line[1];
        X0[2] = X_line[2];
        X0[3] = X_line[3];
        // Only compute the jacobian if there is another Newton iteration to come
        if (norm > CONVERGENCE_NEWTON) {
          double tmp_m[3][3];
          MathExtra::times3_transpose(hessi, Ri, tmp_m);
          MathExtra::times3(Ri, tmp_m, hessi);
          MathExtra::times3_transpose(hessj, Rj, tmp_m);
          MathExtra::times3(Rj, tmp_m, hessj);
          compute_jacobian(gradi, hessi, gradj, hessj, X0[3], jacobian);
        } else {
          converged = true;
          // TODO: consider testing picking the normal with the least error
          //       i.e., likely the grain with the smallest curvature (Hessian norm)
          //       or some other measure like average gradients.
          //       right now we use the gradient on grain i for simplicity and performance. When testing, we could see if using  is just as good
          MathExtra::normalize3(gradi, nij);
        }
        break;
      }
    }

    // If no descent with line search, take full step, try to escape bad region
    if (iter_ls == ITERMAX_LINESEARCH) {
      X0[0] += rhs[0];
      X0[1] += rhs[1];
      X0[2] += rhs[2];
      X0[3] += rhs[3];
      norm = compute_residual_and_jacobian(xci, Ri, shapei, blocki, flagi, xcj, Rj, shapej, blockj, flagj, X0, shapefunc, residual, jacobian);
    }

    if (converged)
      break;
  }

  // LAPACK error are within [-4, 4], use 5 non-touching, -5 non-converging
  if (!converged)
    return -5;
  if (shapefunc[0] > 0.0 || shapefunc[1] > 0.0)
    return 5;

  return 0;
}


int determine_flag(const double* block) {
  const double EPSBLOCK(1e-3);
  int flag(2);
  if ((std::fabs(block[0] - 2) <= EPSBLOCK) && (std::fabs(block[1] - 2) <= EPSBLOCK))
    flag = 0;
  else if (std::fabs(block[0] - block[1]) <= EPSBLOCK)
    flag = 1;
  return flag;
}

// Functions to compute shape function and gradient only when called for newton method
// to avoid computing hessian when not needed and having smoother landscape for the line search
// General case for n1 != n2 > 2
double stable_shape_and_gradient_local_superquad(const double* xlocal, const double* shape, const double* block, double* grad) {
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = std::fabs(xlocal[0] * a_inv);
  double y_b = std::fabs(xlocal[1] * b_inv);
  double z_c = std::fabs(xlocal[2] * c_inv);
  double n1 = block[0];
  double n2 = block[1];
  double x_a_pow_n2_m2 = std::pow(x_a, n2 - 2.0);
  double x_a_pow_n2_m1 = x_a_pow_n2_m2 * x_a;
  double y_b_pow_n2_m2 = std::pow(y_b, n2 - 2.0);
  double y_b_pow_n2_m1 = y_b_pow_n2_m2 * y_b;

  double nu = (x_a_pow_n2_m1 * x_a) + (y_b_pow_n2_m1 * y_b);
  double nu_pow_n1_n2_m2 = std::pow(nu, n1/n2 - 2.0);
  double nu_pow_n1_n2_m1 = nu_pow_n1_n2_m2 * nu;

  double z_c_pow_n1_m2 = std::pow(z_c, n1 -2.0);
  double z_c_pow_n1_m1 = z_c_pow_n1_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n1 * a_inv * x_a_pow_n2_m1 * nu_pow_n1_n2_m1 * signx;
  grad[1] = n1 * b_inv * y_b_pow_n2_m1 * nu_pow_n1_n2_m1 * signy;
  grad[2] = n1 * c_inv * z_c_pow_n1_m1 * signz;

  double F = (nu_pow_n1_n2_m1 * nu) + (z_c_pow_n1_m1 * z_c);

  double scale_factor = std::pow(F, 1.0/n1 -1.0) / n1;

  grad[0] *= scale_factor;
  grad[1] *= scale_factor;
  grad[2] *= scale_factor;

  return (nu_pow_n1_n2_m1 * nu) + (z_c_pow_n1_m1 * z_c) - 1.0;
}

// Special case for n2 = n2 = n > 2
double stable_shape_and_gradient_local_n1equaln2(const double* xlocal, const double* shape, const double n, double* grad) {
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = std::fabs(xlocal[0] * a_inv);
  double y_b = std::fabs(xlocal[1] * b_inv);
  double z_c = std::fabs(xlocal[2] * c_inv);
  double x_a_pow_n_m2 = std::pow(x_a, n - 2.0);
  double x_a_pow_n_m1 = x_a_pow_n_m2 * x_a;
  double y_b_pow_n_m2 = std::pow(y_b, n - 2.0);
  double y_b_pow_n_m1 = y_b_pow_n_m2 * y_b;
  double z_c_pow_n_m2 = std::pow(z_c, n - 2.0);
  double z_c_pow_n_m1 = z_c_pow_n_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n * a_inv * x_a_pow_n_m1 * signx;
  grad[1] = n * b_inv * y_b_pow_n_m1 * signy;
  grad[2] = n * c_inv * z_c_pow_n_m1 * signz;

  double F = (x_a_pow_n_m1 * x_a) + (y_b_pow_n_m1 * y_b) + (z_c_pow_n_m1 * z_c);
  double scale_factor = std::pow(F, 1.0/n -1.0) / n;

  grad[0] *= scale_factor;
  grad[1] *= scale_factor;
  grad[2] *= scale_factor;

  return (x_a_pow_n_m1 * x_a) + (y_b_pow_n_m1 * y_b) + (z_c_pow_n_m1 * z_c) - 1.0;
}


// Special case for n1 = n2 = 2
double stable_shape_and_gradients_local_ellipsoid(const double* xlocal, const double* shape, double* grad) {
  double a = 2.0 / (shape[0] * shape[0]);
  double b = 2.0 / (shape[1] * shape[1]);
  double c = 2.0 / (shape[2] * shape[2]);

  // Equation (14) simplified for n1 = n2 = 2
  grad[0] = a * xlocal[0];
  grad[1] = b * xlocal[1];
  grad[2] = c * xlocal[2];

  return 0.5 * (grad[0]*xlocal[0] + grad[1]*xlocal[1] + grad[2]*xlocal[2]) - 1.0;
}

// Newton Rapson method to find the overlap distance from the contact point given the normal
double compute_overlap_distance(
  const double* shape, const double* block, const double Rot[3][3], const int flag,
  const double* global_point, const double* global_normal,
  const double* center) {
  double local_point[3], local_normal[3];
  double del[3];
  double overlap;
  MathExtra::sub3(global_point, center, del);  // bring origin to 0.0
  MathExtra::transpose_matvec(Rot, del, local_point); 
  MathExtra::transpose_matvec(Rot, global_normal, local_normal);
  
  double local_f;
  double local_grad[3];
  
  // elliposid analytical solution, might need to double check the math 
  // there is an easy way to find this by parametrizing the straight line as
  // X0 + t * n anf then substituting in the ellipsoid equation  for x, y, z
  // this results in a quadratic equation and we take the positive solution since
  // we are taking the outward facing normal for each grain

  if (flag == 0){

    double a_inv2 = 1.0 / (shape[0] * shape[0]);
    double b_inv2 = 1.0 / (shape[1] * shape[1]);
    double c_inv2 = 1.0 / (shape[2] * shape[2]);

    // Coefficients for At^2 + Bt + C = 0
    double A = (local_normal[0] * local_normal[0] * a_inv2) +
               (local_normal[1] * local_normal[1] * b_inv2) +
               (local_normal[2] * local_normal[2] * c_inv2);

    double B = 2.0 * ( (local_point[0] * local_normal[0] * a_inv2) +
                     (local_point[1] * local_normal[1] * b_inv2) +
                     (local_point[2] * local_normal[2] * c_inv2) );

    double C = (local_point[0] * local_point[0] * a_inv2) +
               (local_point[1] * local_point[1] * b_inv2) +
               (local_point[2] * local_point[2] * c_inv2) - 1.0;

    // Discriminant
    double delta = B*B - 4.0*A*C;

    // Clamp delta to zero just in case numerical noise makes it negative
    if (delta < 0.0) delta = 0.0; 
    overlap = (-B + std::sqrt(delta)) / (2.0 * A);
  } else {
      // --- Superquadric Case (Newton-Raphson on Distance Estimator) ---
    
    overlap = 0.0; // Distance along the normal
    double current_p[3];
    double val;
    for (unsigned int iter = 0; iter < ITERMAX_OVERLAP; iter++) {
      // Update current search position: P = Start + t * Normal
      current_p[0] = local_point[0] + overlap * local_normal[0];
      current_p[1] = local_point[1] + overlap * local_normal[1];
      current_p[2] = local_point[2] + overlap * local_normal[2];

      // Calculate Distance Estimator value and Gradient
      if (flag == 1) {
        val = stable_shape_and_gradient_local_n1equaln2(current_p, shape, block[0], local_grad);
      } else {
        val = stable_shape_and_gradient_local_superquad(current_p, shape, block, local_grad);
      }

      // Convergence Check
      if (std::fabs(val) < CONVERGENCE_OVERLAP) break;

      // Newton Step
      double slope = local_grad[0] * local_normal[0] +
                     local_grad[1] * local_normal[1] +
                     local_grad[2] * local_normal[2];

      // Safety check to prevent divide-by-zero if ray grazes surface
      if (std::fabs(slope) < MINSLOPE_OVERLAP) break;

      overlap -= val / slope;
    }
  }
  return overlap;
} 

} // namespace MathExtraSuperellipsoids
