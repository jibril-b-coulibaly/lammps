/* -*- c++ -*- ----------------------------------------------------------
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

#ifndef LMP_MATH_EXTRA_SUPERELLIPOIDS_H
#define LMP_MATH_EXTRA_SUPERELLIPOIDS_H

#include "math_extra.h"

namespace MathExtraSuperellipsoids {
  double beta_func(double a, double b);
  void volume_superellipsoid(const double *blockiness, const double *shape, double volume); // duplicated from math_extra might remove
  void inertia_superellipsoid(const double *shape, const double *blockiness, double density, double *inertia); // duplicated from math_extra might remove

  // needed for shape functions grad and matrix 
  void local2global_vector(const double v[3], const double *quat, double global_v[3]);
  void global2local_vector(const double v[3], const double *quat, double local_v[3]);
  void local2global_matrix(const double m[3][3], const double *quat, double global_m[3][3]);
  void global2local_matrix(const double m[3][3], const double *quat, double local_m[3][3]);

  // shape function computations
  void shape_function_local(const double *shape, const double *block, const double *quat, const double *point, double local_f);
  void shape_function_global(const double *shape, const double *block, const double *quat, const double *point, double global_f);
  void shape_function_local_grad(const double *shape, const double *block, const double *quat, const double *point, double *local_grad);
  void shape_function_local_hessian(const double *shape, const double *block, const double *quat, const double *point, double local_hessian[3][3]);

  inline double det4_M44_zero(const double m[4][4]);


  // ADD CONTACT DETECTION HERE

};


/* ----------------------------------------------------------------------
   determinant of a 4x4 matrix M with M[3][3] assumed to be zero
------------------------------------------------------------------------- */
inline double MathExtraSuperellipsoids::det4_M44_zero(const double m[4][4])
{
    // Define the 3x3 submatrices (M_41, M_42, M_43)

    // Submatrix M_41 
    double m41[3][3] = {
        {m[0][1], m[0][2], m[0][3]},
        {m[1][1], m[1][2], m[1][3]},
        {m[2][1], m[2][2], m[2][3]}
    };

    // Submatrix M_42 
    double m42[3][3] = {
        {m[0][0], m[0][2], m[0][3]},
        {m[1][0], m[1][2], m[1][3]},
        {m[2][0], m[2][2], m[2][3]}
    };

    // Submatrix M_43
    double m43[3][3] = {
        {m[0][0], m[0][1], m[0][3]},
        {m[1][0], m[1][1], m[1][3]},
        {m[2][0], m[2][1], m[2][3]}
    };
    
    // Calculate the determinant using the simplified Laplace expansion (M_44=0)
    // det(M) = -M[3][0]*det(M_41) + M[3][1]*det(M_42) - M[3][2]*det(M_43)
    
    double ans = -m[3][0] * MathExtra::det3(m41) 
                 + m[3][1] * MathExtra::det3(m42) 
                 - m[3][2] * MathExtra::det3(m43);
                 
    return ans;
}




#endif
