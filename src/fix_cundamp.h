/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Damping used in Yade-DEM, reduction of unbalanced force
   written by Jibril B. Coulibaly @ Northwestern University, 04/12/2019
   Adapted for rigid bodies: damping applied to the rigid body, not individual particles
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(cundamp,FixCundamp)

#else

#ifndef LMP_FIX_CUNDAMP_H
#define LMP_FIX_CUNDAMP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCundamp : public Fix {
 public:
  FixCundamp(class LAMMPS *, int, char **);
  virtual ~FixCundamp();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

 protected:
  double *gamma_linear,*gamma_angular;
  int ilevel_respa;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
