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

namespace MathExtraSuperellipsoids {

/* ----------------------------------------------------------------------
   beta function B(x,y) = Gamma(x) * Gamma(y) / Gamma(x+y)
------------------------------------------------------------------------- */
double beta_func(double a, double b) {
    return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

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


}


