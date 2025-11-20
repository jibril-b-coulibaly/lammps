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

#include <cmath>
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

  // 4 by 4 sytems solvers, they all overwrite b with the solution
  inline bool solve_4x4_manual(double A[16], double b[4]);  
  inline bool solve_4x4_robust(double A[16], double b[4]);
  inline bool solve_4x4_robust_unrolled(double A[16], double b[4]); 

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

inline bool MathExtraSuperellipsoids::solve_4x4_manual(double A[16], double b[4]) {
    // 1. Pivot 0 
    double inv0 = 1.0 / A[0];
    double m1 = A[4] * inv0;
    double m2 = A[8] * inv0;
    double m3 = A[12] * inv0;

    A[5] -= m1 * A[1]; A[6] -= m1 * A[2]; A[7] -= m1 * A[3]; b[1] -= m1 * b[0];
    A[9] -= m2 * A[1]; A[10] -= m2 * A[2]; A[11] -= m2 * A[3]; b[2] -= m2 * b[0];
    A[13] -= m3 * A[1]; A[14] -= m3 * A[2]; A[15] -= m3 * A[3]; b[3] -= m3 * b[0];

    // 2. Pivot 1 
    double inv1 = 1.0 / A[5];
    double m4 = A[9] * inv1;
    double m5 = A[13] * inv1;

    A[10] -= m4 * A[6]; A[11] -= m4 * A[7]; b[2] -= m4 * b[1];
    A[14] -= m5 * A[6]; A[15] -= m5 * A[7]; b[3] -= m5 * b[1];

    // 3. Pivot 2
    double inv2 = 1.0 / A[10];
    double m6 = A[14] * inv2;

    A[15] -= m6 * A[11]; b[3] -= m6 * b[2];

    // 4. Backward Substitution
    b[3] = b[3] / A[15];
    b[2] = (b[2] - A[11] * b[3]) * inv2;
    b[1] = (b[1] - A[7] * b[3] - A[6] * b[2]) * inv1;
    b[0] = (b[0] - A[3] * b[3] - A[2] * b[2] - A[1] * b[1]) * inv0;

    return true;
}

inline bool MathExtraSuperellipsoids::solve_4x4_robust(double A[16], double b[4]) {
    // Helper lambda to access A[row, col]
    auto at = [&](int r, int c) -> double& { return A[r * 4 + c]; };

    // --- FORWARD ELIMINATION with PARTIAL PIVOTING ---
    
    for (int i = 0; i < 3; ++i) { // Loop over columns 0, 1, 2
        // 1. Find the Pivot (Max absolute value in this column)
        int pivot_row = i;
        double max_val = std::abs(at(i, i));

        for (int k = i + 1; k < 4; ++k) {
            double val = std::abs(at(k, i));
            if (val > max_val) {
                max_val = val;
                pivot_row = k;
            }
        }

        // 2. Singularity Check (The "Flat Particle" Guard)
        if (max_val < 1e-14) return false;

        // 3. Swap Rows if needed (Swap A rows AND b elements)
        if (pivot_row != i) {
            std::swap(b[i], b[pivot_row]);
            for (int k = i; k < 4; ++k) { // Only need to swap from column 'i' onwards
                std::swap(at(i, k), at(pivot_row, k));
            }
        }

        // 4. Eliminate
        double inv_pivot = 1.0 / at(i, i);
        for (int k = i + 1; k < 4; ++k) {
            double factor = at(k, i) * inv_pivot;
            // A[k, i] becomes 0, no need to compute it.
            // Update the rest of the row:
            for (int j = i + 1; j < 4; ++j) {
                at(k, j) -= factor * at(i, j);
            }
            // Update RHS
            b[k] -= factor * b[i];
        }
    }

    // Final Pivot Check for the last element
    if (std::abs(at(3, 3)) < 1e-14) return false;

    // --- BACKWARD SUBSTITUTION ---
    b[3] /= at(3, 3);
    b[2] = (b[2] - at(2, 3) * b[3]) / at(2, 2);
    b[1] = (b[1] - at(1, 2) * b[2] - at(1, 3) * b[3]) / at(1, 1);
    b[0] = (b[0] - at(0, 1) * b[1] - at(0, 2) * b[2] - at(0, 3) * b[3]) / at(0, 0);

    return true;

}

inline bool MathExtraSuperellipsoids::solve_4x4_robust_unrolled(double A[16], double b[4]) {
     // --- COLUMN 0 ---
    // 1. Find Pivot in Col 0
    int p = 0; 
    double max_val = std::abs(A[0]);
    double val;

    val = std::abs(A[4]); 
    if (val > max_val) { max_val = val; p = 1; }
    val = std::abs(A[8]); 
    if (val > max_val) { max_val = val; p = 2; }
    val = std::abs(A[12]);
    if (val > max_val) { max_val = val; p = 3; }

    if (max_val < 1e-14) return false;
    // 2. Swap Row 0 with Row p
    if (p != 0) {
        int row_offset = p * 4;
        std::swap(b[0], b[p]);
        std::swap(A[0], A[row_offset]);     std::swap(A[1], A[row_offset + 1]);
        std::swap(A[2], A[row_offset + 2]); std::swap(A[3], A[row_offset + 3]);
    }

    // 3. Eliminate Col 0
    {
        double inv = 1.0 / A[0];
        // Row 1
        double f1 = A[4] * inv;
        A[5] -= f1 * A[1]; A[6] -= f1 * A[2]; A[7] -= f1 * A[3]; b[1] -= f1 * b[0];
        // Row 2
        double f2 = A[8] * inv;
        A[9] -= f2 * A[1]; A[10] -= f2 * A[2]; A[11] -= f2 * A[3]; b[2] -= f2 * b[0];
        // Row 3
        double f3 = A[12] * inv;
        A[13] -= f3 * A[1]; A[14] -= f3 * A[2]; A[15] -= f3 * A[3]; b[3] -= f3 * b[0];
    }

    // --- COLUMN 1 ---
    // 1. Find Pivot in Col 1 (starting from row 1)
    p = 1;
    max_val = std::abs(A[5]);
    
    val = std::abs(A[9]);  if (val > max_val) { max_val = val; p = 2; }
    val = std::abs(A[13]); if (val > max_val) { max_val = val; p = 3; }

    if (max_val < 1e-14) return false;

    // 2. Swap Row 1 with Row p
    if (p != 1) {
        int row_offset = p * 4;
        std::swap(b[1], b[p]);
        // Optimization: Col 0 is already 0, so we only swap cols 1,2,3
        std::swap(A[5], A[row_offset + 1]);
        std::swap(A[6], A[row_offset + 2]);
        std::swap(A[7], A[row_offset + 3]);
    }

    // 3. Eliminate Col 1
    {
        double inv = 1.0 / A[5];
        // Row 2
        double f2 = A[9] * inv;
        A[10] -= f2 * A[6]; A[11] -= f2 * A[7]; b[2] -= f2 * b[1];
        // Row 3
        double f3 = A[13] * inv;
        A[14] -= f3 * A[6]; A[15] -= f3 * A[7]; b[3] -= f3 * b[1];
    }

    // --- COLUMN 2 ---
    // 1. Find Pivot in Col 2 (starting from row 2)
    p = 2;
    max_val = std::abs(A[10]);

    val = std::abs(A[14]); if (val > max_val) { max_val = val; p = 3; }

    if (max_val < 1e-14) return false;

    // 2. Swap Row 2 with Row p
    if (p != 2) {
        std::swap(b[2], b[3]);
        // Optimization: Only swap cols 2,3
        std::swap(A[10], A[14]);
        std::swap(A[11], A[15]);
    }

    // 3. Eliminate Col 2
    {
        double inv = 1.0 / A[10];
        // Row 3
        double f3 = A[14] * inv;
        A[15] -= f3 * A[11]; b[3] -= f3 * b[2];
    }

    // --- BACKWARD SUBSTITUTION ---
    // Check last pivot
    if (std::abs(A[15]) < 1e-14) return false;

    double inv3 = 1.0 / A[15];
    b[3] *= inv3;

    double inv2 = 1.0 / A[10];
    b[2] = (b[2] - A[11] * b[3]) * inv2;

    double inv1 = 1.0 / A[5];
    b[1] = (b[1] - A[6] * b[2] - A[7] * b[3]) * inv1;

    double inv0 = 1.0 / A[0];
    b[0] = (b[0] - A[1] * b[1] - A[2] * b[2] - A[3] * b[3]) * inv0;

    return true;

}



#endif
