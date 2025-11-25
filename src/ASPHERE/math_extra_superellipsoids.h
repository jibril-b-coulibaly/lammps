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
  inline constexpr double TIKHONOV_SCALE = 1e-8;
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
  inline bool check_oriented_bounding_boxes(const double* xc1, const double R1[3][3], const double* shape1,
                                        const double* xc2, const double R2[3][3], const double* shape2, 
                                        double* cached_axis);

  inline bool check_intersection_axis(const int axis_id, const double C[3][3], const double AbsC[3][3], 
                                      const double* center_distance_box1, const double* center_distance_box2,
                                      const double* a, const double* b);

  inline bool check_collision_and_get_seed(const double* xc1, const double R1[3][3], const double* shape1,
                                         const double* xc2, const double R2[3][3], const double* shape2,
                                        double* cached_axis, double* contact_point);


  // Jibril's versions of the functions for contact detection
  double shape_and_derivatives_local(const double* xlocal, const double* shape, const double* block, const int flag, double* grad, double hess[3][3]);
  double shape_and_derivatives_local_superquad(const double* xlocal, const double* shape, const double* block, double* grad, double hess[3][3]);
  double shape_and_derivatives_local_n1equaln2(const double* xlocal, const double* shape, const double n, double* grad, double hess[3][3]);
  double shape_and_derivatives_local_ellipsoid(const double* xlocal, const double* shape, double* grad, double hess[3][3]);
  double shape_and_derivatives_global(const double* xc, const double R[3][3], const double* shape, const double* block, const int flag, const double* X0, double* grad, double hess[3][3]);
  double compute_residual(const double shapefunci, const double* gradi_global, const double shapefuncj, const double* gradj_global, const double mu2, double* residual);
  void compute_jacobian(const double* gradi_global, const double hessi_global[3][3], const double* gradj_global, const double hessj_global[3][3], const double mu2, double* jacobian);
  double compute_residual_and_jacobian(const double* xci, const double Ri[3][3], const double* shapei, const double* blocki, const int flagi,
                                       const double* xcj, const double Rj[3][3], const double* shapej, const double* blockj, const int flagj,
                                       const double* X, double* shapefunc, double* residual, double* jacobian);
  int determine_contact_point(const double* xci, const double Ri[3][3], const double* shapei, const double* blocki, const int flagi,
                              const double* xcj, const double Rj[3][3], const double* shapej, const double* blockj, const int flagj,
                              double* X0, double* nij);
  int determine_flag(const double* block);
 
  // functions to compute shape function and gradient only when called for newton method
  double stable_shape_and_gradient_local_superquad(const double* xlocal, const double* shape, const double* block, double* grad);
  double stable_shape_and_gradient_local_n1equaln2(const double* xlocal, const double* shape, const double n, double* grad);
  double stable_shape_and_gradient_local_ellipsoid(const double* xlocal, const double* shape, double* grad);
  double compute_overlap_distance(const double* shape, const double* block, const double Rot[3][3], const int flag, const double* global_point, const double* global_normal, const double* center);

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
    
    // Tikhonov regularization
    // High blockiness grains can have zero curvature / singular Hessian
    // along principal local axes (x=0, y=0, z=0)
    const double diag_weight = TIKHONOV_SCALE * (A[0] + A[5] + A[10] + A[15]);
    A[0]  += diag_weight;
    A[5]  += diag_weight;
    A[10] += diag_weight;
    A[15] += diag_weight;

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

    // Tikhonov regularization
    // High blockiness grains can have zero curvature / singular Hessian
    // along principal local axes (x=0, y=0, z=0)
    const double diag_weight = TIKHONOV_SCALE * (A[0] + A[5] + A[10] + A[15]);
    A[0]  += diag_weight;
    A[5]  += diag_weight;
    A[10] += diag_weight;
    A[15] += diag_weight;

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
    
    // Tikhonov regularization
    // High blockiness grains can have zero curvature / singular Hessian
    // along principal local axes (x=0, y=0, z=0)
    const double diag_weight = TIKHONOV_SCALE * (A[0] + A[5] + A[10] + A[15]);
    A[0]  += diag_weight;
    A[5]  += diag_weight;
    A[10] += diag_weight;
    A[15] += diag_weight;
    
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


// algorithm from https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
/* * Oriented Bounding Box intersection test.
 * Logic and optimization strategies adapted from LIGGGHTS (CFDEMproject).
 * See: src/math_extra_liggghts_nonspherical.cpp in LIGGGHTS distribution.
 * * This implementation uses the "cached separating axis" optimization 
 * for temporal coherence.
 */
inline bool MathExtraSuperellipsoids::check_oriented_bounding_boxes(
    const double* xc1, const double R1[3][3], const double* shape1,
    const double* xc2, const double R2[3][3], const double* shape2, 
    double* cached_axis
){
    // cache axis is the axis that separated the boxes last time
    // due to temporal coherence we check it first

    bool separated = false;

    // for orientated bounding boxes we check the 15 separating axes
    double C[3][3], AbsC[3][3];
    MathExtra::transpose_times3(R1, R2, C); // C = R1^T * R2
    for (unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<3; j++){
            AbsC[i][j] = std::fabs(C[i][j]); // for when absolute values are needed
        }
    }

    double center_distance[3];
    for (unsigned int i=0; i<3; i++){
        center_distance[i] = xc2[i] - xc1[i];
    } 

    // Project center distance into both local frames
    double center_distance_box1[3], center_distance_box2[3];
    MathExtra::transpose_matvec(R1, center_distance,  center_distance_box1);
    MathExtra::transpose_matvec(R2, center_distance,  center_distance_box2);

    // first check the cached axis
    const int axis = (int) (*cached_axis);
    separated = check_intersection_axis(axis, C, AbsC, center_distance_box1, center_distance_box2, shape1, shape2);

    if (separated) return true;
    // then check all the other axes
    for (int axis_id = 0; axis_id < 15; axis_id++){
        if (axis_id == axis) continue; // already checked
        separated = check_intersection_axis(axis_id, C, AbsC, center_distance_box1, center_distance_box2, shape1, shape2);
        if (separated) {
            *cached_axis = axis_id; // update cached axis
            return true;
        }
    }
    return false; // no separation found
}

inline bool MathExtraSuperellipsoids::check_intersection_axis(
    const int axis_id, const double C[3][3], const double AbsC[3][3], 
    const double* center_distance_box1, const double* center_distance_box2,
    const double* a, const double* b
){
    // here axis_id goes from 0 to 14
    // a and b are the half-sizes of the boxes along their local axes
    // returns true if there is a separation along this axis
    // changes the cached axis if separation found
    double R1, R2, R;

    switch(axis_id){
        case 0: // A0
            R1 = a[0];
            R2 = b[0] * AbsC[0][0] + b[1] * AbsC[0][1] + b[2] * AbsC[0][2];
            R = std::fabs(center_distance_box1[0]);
            break;
        case 1: // A1
            R1 = a[1];
            R2 = b[0] * AbsC[1][0] + b[1] * AbsC[1][1] + b[2] * AbsC[1][2];
            R = std::fabs(center_distance_box1[1]);
            break;
        case 2: // A2
            R1 = a[2];
            R2 = b[0] * AbsC[2][0] + b[1] * AbsC[2][1] + b[2] * AbsC[2][2];
            R = std::fabs(center_distance_box1[2]);
            break;
        case 3: // B0
            R1 = a[0] * AbsC[0][0] + a[1] * AbsC[1][0] + a[2] * AbsC[2][0];
            R2 = b[0];
            R = std::fabs(center_distance_box2[0]);
            break;
        case 4: // B1
            R1 = a[0] * AbsC[0][1] + a[1] * AbsC[1][1] + a[2] * AbsC[2][1];
            R2 = b[1];
            R = std::fabs(center_distance_box2[1]);
            break;
        case 5: // B2
            R1 = a[0] * AbsC[0][2] + a[1] * AbsC[1][2] + a[2] * AbsC[2][2];
            R2 = b[2];
            R = std::fabs(center_distance_box2[2]);
            break;
        case 6: // A0 x B0
            R1 = a[1] * AbsC[2][0] + a[2] * AbsC[1][0];
            R2 = b[1] * AbsC[0][2] + b[2] * AbsC[0][1];
            R = std::fabs(center_distance_box1[2] * C[1][0] - center_distance_box1[1] * C[2][0]);
            break;
        case 7: // A0 x B1
            R1 = a[1] * AbsC[2][1] + a[2] * AbsC[1][1];
            R2 = b[0] * AbsC[0][2] + b[2] * AbsC[0][0];
            R = std::fabs(center_distance_box1[2] * C[1][1] - center_distance_box1[1] * C[2][1]);
            break;
        case 8: // A0 x B2
            R1 = a[1] * AbsC[2][2] + a[2] * AbsC[1][2];
            R2 = b[0] * AbsC[0][1] + b[1] * AbsC[0][0];
            R = std::fabs(center_distance_box1[2] * C[1][2] - center_distance_box1[1] * C[2][2]);
            break;
        case 9: // A1 x B0
            R1 = a[0] * AbsC[2][0] + a[2] * AbsC[0][0];
            R2 = b[1] * AbsC[1][2] + b[2] * AbsC[1][1];
            R = std::fabs(center_distance_box1[0] * C[2][0] - center_distance_box1[2] * C[0][0]);
            break;
        case 10: // A1 x B1
            R1 = a[0] * AbsC[2][1] + a[2] * AbsC[0][1];
            R2 = b[0] * AbsC[1][2] + b[2] * AbsC[1][0];
            R = std::fabs(center_distance_box1[0] * C[2][1] - center_distance_box1[2] * C[0][1]);
            break;
        case 11: // A1 x B2
            R1 = a[0] * AbsC[2][2] + a[2] * AbsC[0][2];
            R2 = b[0] * AbsC[1][1] + b[1] * AbsC[1][0];
            R = std::fabs(center_distance_box1[0] * C[2][2] - center_distance_box1[2] * C[0][2]);
            break;
        case 12: // A2 x B0
            R1 = a[0] * AbsC[1][0] + a[1] * AbsC[0][0];
            R2 = b[1] * AbsC[2][2] + b[2] * AbsC[2][1];
            R = std::fabs(center_distance_box1[1] * C[0][0] - center_distance_box1[0] * C[1][0]);
            break;
        case 13: // A2 x B1
            R1 = a[0] * AbsC[1][1] + a[1] * AbsC[0][1];
            R2 = b[0] * AbsC[2][2] + b[2] * AbsC[2][0];
            R = std::fabs(center_distance_box1[1] * C[0][1] - center_distance_box1[0] * C[1][1]);
            break;
        case 14: // A2 x B2
            R1 = a[0] * AbsC[1][2] + a[1] * AbsC[0][2];
            R2 = b[0] * AbsC[2][1] + b[1] * AbsC[2][0];
            R = std::fabs(center_distance_box1[1] * C[0][2] - center_distance_box1[0] * C[1][2]);
            break;
    }

    if (R > R1 + R2){
        return true; // separation found
    } else {
        return false; // no separation
    }
}


inline bool MathExtraSuperellipsoids::check_collision_and_get_seed(
    const double* xc1, const double R1[3][3], const double* shape1,
    const double* xc2, const double R2[3][3], const double* shape2, 
    double* cached_axis, double* contact_point
){  
    // cache axis is the axis that separated the boxes last time
    // due to temporal coherence we check it first

    double C[3][3], AbsC[3][3];
    MathExtra::transpose_times3(R1, R2, C); // C = R1^T * R2
    
    // for orientated bounding boxes we check the 15 separating axes
    const double eps = 1e-20;
    for (unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<3; j++){
            AbsC[i][j] = std::fabs(C[i][j]) + eps; // Add epsilon to prevent division by zero in edge cases
        }
    }

    double center_distance[3]; // Center distance in Global Frame
        for (unsigned int i=0; i<3; i++){
        center_distance[i] = xc2[i] - xc1[i];
    } 

    // Project center distance into both local frames
    double center_distance_box1[3], center_distance_box2[3];
    MathExtra::transpose_matvec(R1, center_distance, center_distance_box1);
    MathExtra::transpose_matvec(R2, center_distance, center_distance_box2);

    int best_axis = -1;
    double min_overlap = 0.0;
    const double edge_bias = 1.05; // Prefer face contacts over edge contacts

    // Lambda to test an axis. Returns TRUE if SEPARATED.
    // I was reading that lambdas can be optimized away by the compiler.
    // and have less overhead than function calls.
    auto test_axis_separated = [&](int i) -> bool {
        double R1_rad, R2_rad, dist, overlap;

        // Switch is efficient here; compiler generates a jump table.
        switch(i){
            case 0: // A0
                R1_rad = shape1[0];
                R2_rad = shape2[0] * AbsC[0][0] + shape2[1] * AbsC[0][1] + shape2[2] * AbsC[0][2];
                dist = std::fabs(center_distance_box1[0]);
                break;
            case 1: // A1
                R1_rad = shape1[1];
                R2_rad = shape2[0] * AbsC[1][0] + shape2[1] * AbsC[1][1] + shape2[2] * AbsC[1][2];
                dist = std::fabs(center_distance_box1[1]);
                break;
            case 2: // A2
                R1_rad = shape1[2];
                R2_rad = shape2[0] * AbsC[2][0] + shape2[1] * AbsC[2][1] + shape2[2] * AbsC[2][2];
                dist = std::fabs(center_distance_box1[2]);
                break;
            case 3: // B0
                R1_rad = shape1[0] * AbsC[0][0] + shape1[1] * AbsC[1][0] + shape1[2] * AbsC[2][0];
                R2_rad = shape2[0];
                dist = std::fabs(center_distance_box2[0]);
                break;
            case 4: // B1
                R1_rad = shape1[0] * AbsC[0][1] + shape1[1] * AbsC[1][1] + shape1[2] * AbsC[2][1];
                R2_rad = shape2[1];
                dist = std::fabs(center_distance_box2[1]);
                break;
            case 5: // B2
                R1_rad = shape1[0] * AbsC[0][2] + shape1[1] * AbsC[1][2] + shape1[2] * AbsC[2][2];
                R2_rad = shape2[2];
                dist = std::fabs(center_distance_box2[2]);
                break;
            case 6: // A0 x B0
                R1_rad = shape1[1] * AbsC[2][0] + shape1[2] * AbsC[1][0];
                R2_rad = shape2[1] * AbsC[0][2] + shape2[2] * AbsC[0][1];
                dist = std::fabs(center_distance_box1[2] * C[1][0] - center_distance_box1[1] * C[2][0]);
                break;
            case 7: // A0 x B1
                R1_rad = shape1[1] * AbsC[2][1] + shape1[2] * AbsC[1][1];
                R2_rad = shape2[0] * AbsC[0][2] + shape2[2] * AbsC[0][0];
                dist = std::fabs(center_distance_box1[2] * C[1][1] - center_distance_box1[1] * C[2][1]);
                break;
            case 8: // A0 x B2
                R1_rad = shape1[1] * AbsC[2][2] + shape1[2] * AbsC[1][2];
                R2_rad = shape2[0] * AbsC[0][1] + shape2[1] * AbsC[0][0];
                dist = std::fabs(center_distance_box1[2] * C[1][2] - center_distance_box1[1] * C[2][2]);
                break;
            case 9: // A1 x B0
                R1_rad = shape1[0] * AbsC[2][0] + shape1[2] * AbsC[0][0];
                R2_rad = shape2[1] * AbsC[1][2] + shape2[2] * AbsC[1][1];
                dist = std::fabs(center_distance_box1[0] * C[2][0] - center_distance_box1[2] * C[0][0]);
                break;
            case 10: // A1 x B1
                R1_rad = shape1[0] * AbsC[2][1] + shape1[2] * AbsC[0][1];
                R2_rad = shape2[0] * AbsC[1][2] + shape2[2] * AbsC[1][0];
                dist = std::fabs(center_distance_box1[0] * C[2][1] - center_distance_box1[2] * C[0][1]);
                break;
            case 11: // A1 x B2
                R1_rad = shape1[0] * AbsC[2][2] + shape1[2] * AbsC[0][2];
                R2_rad = shape2[0] * AbsC[1][1] + shape2[1] * AbsC[1][0];
                dist = std::fabs(center_distance_box1[0] * C[2][2] - center_distance_box1[2] * C[0][2]);
                break;
            case 12: // A2 x B0
                R1_rad = shape1[0] * AbsC[1][0] + shape1[1] * AbsC[0][0];
                R2_rad = shape2[1] * AbsC[2][2] + shape2[2] * AbsC[2][1];
                dist = std::fabs(center_distance_box1[1] * C[0][0] - center_distance_box1[0] * C[1][0]);
                break;
            case 13: // A2 x B1
                R1_rad = shape1[0] * AbsC[1][1] + shape1[1] * AbsC[0][1];
                R2_rad = shape2[0] * AbsC[2][2] + shape2[2] * AbsC[2][0];
                dist = std::fabs(center_distance_box1[1] * C[0][1] - center_distance_box1[0] * C[1][1]);
                break;
            case 14: // A2 x B2
                R1_rad = shape1[0] * AbsC[1][2] + shape1[1] * AbsC[0][2];
                R2_rad = shape2[0] * AbsC[2][1] + shape2[1] * AbsC[2][0];
                dist = std::fabs(center_distance_box1[1] * C[0][2] - center_distance_box1[0] * C[1][2]);
                break;
            default: return false;
        }

        if (dist > R1_rad + R2_rad) return true; // Separated!

        // If not separated, track the overlap depth
        overlap = (R1_rad + R2_rad) - dist;
        
        // Bias: Penalize edge axes slightly to prefer stable face contacts
        if (i >= 6) overlap *= edge_bias;

        if (overlap < min_overlap) {
            min_overlap = overlap;
            best_axis = i;
        }
        return false; // Not separated
    };

    // Check Cached Axis First (Temporal Coherence)
    int c_axis = (int)(*cached_axis);
    if (test_axis_separated(c_axis)) return false; 

    // Check remaining axes
    for (int i = 0; i < 15; i++){
        if (i == c_axis) continue;
        if (test_axis_separated(i)) {
            *cached_axis = (double)i;
            return false;
        }
    }
   
    // If we reached here, 'best_axis' holds the axis index where the overlap is minimal
    if (best_axis < 6) {
        // Face-to-Face contact logic: Project "Incident" box onto "Reference" face, clip to find overlap center.
        // Pointers to define who is Reference (the face) and who is Incident
        const double* posRef = xc1;
        const double* posInc = xc2;
        const double (*RRef)[3] = R1;
        const double (*RInc)[3] = R2;
        const double* shapeRef = shape1;
        const double* shapeInc = shape2;
        double* D_local_Ref = center_distance_box1; // Center dist in Ref frame

        int axis = best_axis; 

        // Swap if Reference is Box 2 (Indices 3, 4, 5)
        if (best_axis >= 3) {
            posRef = xc2;
            posInc = xc1;
            RRef = R2;
            RInc = R1;
            shapeRef = shape2;
            shapeInc = shape1;
            D_local_Ref = center_distance_box2;
            axis -= 3;
        }

        double seed_local[3];

        //Normal Component: Midway through the penetration depth
        // Calculate projected radius of Incident block onto this axis
        
        double dir = (D_local_Ref[axis] > 0) ? 1.0 : -1.0;
        double radInc_proj = 0.0;
        for(int k=0; k<3; k++) {
            // If swapped (Box 2 is Ref), we need AbsC^T, so we swap AbsC indices
            double val = (best_axis < 3) ? AbsC[axis][k] : AbsC[k][axis];
            radInc_proj += shapeInc[k] * val;
        }

        double surfRef = dir * shapeRef[axis];
        double surfInc = D_local_Ref[axis] - (dir * radInc_proj);
        seed_local[axis] = 0.5 * (surfRef + surfInc);

        // Lateral Components: 1D Interval Overlap
        for(int k=0; k<3; k++) {
            if (k == axis) continue; // Skip the normal axis

            double minRef = -shapeRef[k];
            double maxRef =  shapeRef[k];

            double radInc = 0.0;
            for(int j=0; j<3; j++) {
                double val = (best_axis < 3) ? AbsC[k][j] : AbsC[j][k]; 
                radInc += shapeInc[j] * val;
            }
            double centerInc = D_local_Ref[k];
            
            double minInc = centerInc - radInc;
            double maxInc = centerInc + radInc;

            // Find intersection of intervals [minRef, maxRef] and [minInc, maxInc]
            double start = (minRef > minInc) ? minRef : minInc; 
            double end   = (maxRef < maxInc) ? maxRef : maxInc; 
            seed_local[k] = 0.5 * (start + end); // Midpoint of overlap

        }

        // Transform Local Seed -> World Space
        MathExtra::matvec(RRef, seed_local, contact_point);
        for(int k=0; k<3; k++) contact_point[k] += posRef[k];
    } 
    else {
        // Edge-to-edge contact logic: Midpoint of the closest points on the two skew edge lines.
        // The logic is that index 6 corresponds to A_0 x B_0, 7 to A_0 x B_1, ..., 14 to A_2 x B_2
        int edgeA_idx = (best_axis - 6) / 3;
        int edgeB_idx = (best_axis - 6) % 3;

        // Get World directions of the edges
        double u[3] = { R1[0][edgeA_idx], R1[1][edgeA_idx], R1[2][edgeA_idx] };
        double v[3] = { R2[0][edgeB_idx], R2[1][edgeB_idx], R2[2][edgeB_idx] };

        // Identify the specific edges by checking the normal direction
        // The normal N is roughly the distance vector center_distance for the closest edges
        double N_loc1[3], N_loc2[3];
        MathExtra::transpose_matvec(R1, center_distance, N_loc1);
        MathExtra::transpose_matvec(R2, center_distance, N_loc2);

        // Find Center of Edge A in World Space
        double midA[3]; for(int k=0; k<3; k++) midA[k] = xc1[k];
        for(int k=0; k<3; k++){
            if(k == edgeA_idx) continue;
            // Move to the face pointing towards B
            double sign = (N_loc1[k] > 0) ? 1.0 : -1.0;
            double offset = sign * shape1[k];
            midA[0] += R1[0][k]*offset; midA[1] += R1[1][k]*offset; midA[2] += R1[2][k]*offset;
        }

        // Find Center of Edge B in World Space
        double midB[3]; for(int k=0; k<3; k++) midB[k] = xc2[k];
        for(int k=0; k<3; k++){
            if(k == edgeB_idx) continue;
            // Move to the face pointing away from A (Since center_distance is A->B, we check -N_loc2)
            double sign = (N_loc2[k] < 0) ? 1.0 : -1.0; 
            double offset = sign * shape2[k];
            midB[0] += R2[0][k]*offset; midB[1] += R2[1][k]*offset; midB[2] += R2[2][k]*offset;
        }

        // Closest Points on Two Skew Lines 
        // Line1 parameterized by s: P_A = midA + s*u
        // Line2 parameterized by t: P_B = midB + t*v
        double r[3] = { midB[0]-midA[0], midB[1]-midA[1], midB[2]-midA[2] };
        double u_dot_v = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
        double u_dot_r = u[0]*r[0]+u[1]*r[1]+u[2]*r[2];
        double v_dot_r = v[0]*r[0]+v[1]*r[1]+v[2]*r[2];
        
        // Denom is 1 - (u.v)^2 because u and v are unit vectors
        double denom = 1.0 - u_dot_v*u_dot_v + eps; 
        double s = (u_dot_r - u_dot_v * v_dot_r) / denom;
        double t = (u_dot_v * u_dot_r - v_dot_r) / denom; // Note: simplified derivation

        // Compute World Points
        double PA[3] = { midA[0]+s*u[0], midA[1]+s*u[1], midA[2]+s*u[2] };
        double PB[3] = { midB[0]+t*v[0], midB[1]+t*v[1], midB[2]+t*v[2] };

        // Seed is the midpoint
        for(int k=0; k<3; k++) contact_point[k] = 0.5 * (PA[k] + PB[k]);
    }

    return true; // Collision confirmed
}


#endif
