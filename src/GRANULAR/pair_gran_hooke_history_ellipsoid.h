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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(gran/hooke/history/ellipsoid,PairGranHookeHistoryEllipsoid);
// clang-format on
#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_ELLIPSOID_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_ELLIPSOID_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGranHookeHistoryEllipsoid : public Pair {
 public:
  PairGranHookeHistoryEllipsoid(class LAMMPS *);
  ~PairGranHookeHistoryEllipsoid() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void reset_dt() override;
  double single(int, int, int, int, double, double, double, double &) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  double atom2cut(int) override;
  double radii2cut(double, double) override;

 protected:
  double kn, kt, gamman, gammat, xmu;
  int dampflag;
  double dt;
  int freeze_group_bit;
  int history;
  int limit_damping;

  int neighprev;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;

  int size_history;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, null pointer if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  void allocate();

 private:
  static double shape_and_derivatives_local(const double*, const double*, const double*, const int, double*, double[3][3]);
  static double shape_and_derivatives_local_superquad(const double*, const double*, const double*, double*, double[3][3]);
  static double shape_and_derivatives_local_n1equaln2(const double*, const double*, const double, double*, double[3][3]);
  static double shape_and_derivatives_local_ellipsoid(const double*, const double*, double*, double[3][3]);
  static double shape_and_derivatives_global(const double*, const double[3][3], const double*, const double*, const int, const double*, double*, double[3][3]);

  static double compute_residual(const double, const double*, const double, const double*, const double, double*);
  static void compute_jacobian(const double*, const double[3][3], const double*, const double[3][3], const double, double* jacobian);
  static double compute_residual_and_jacobian(const double*, const double[3][3], const double*, const double*, const int, const double*, const double[3][3], const double*, const double*, const int, const double*, double*, double*, double*);
  static int determine_contact_point(const double*, const double[3][3], const double*, const double*, const double*, const double[3][3], const double*, const double*, double*);
  static int determine_flag(const double*);

  // Below not implemented. Placeholder if we decide not to compute local hessian in line search
  static double shape_and_gradient_local(const double*, const double*, const double*, double*); // would return a vector of temporary variables
  static double hessian_local(const double*, const double*, const double*, double*); // would use the above vector of temporary variables to compute local hessian
};

}    // namespace LAMMPS_NS

#endif
#endif
