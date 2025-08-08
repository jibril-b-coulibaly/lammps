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
PairStyle(gran/hooke/history/ellipsoid,PairGranHookeHistory);
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
};

}    // namespace LAMMPS_NS

#endif
#endif
