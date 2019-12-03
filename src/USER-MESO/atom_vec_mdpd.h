/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(mdpd,AtomVecMDPD)

#else

#ifndef LMP_ATOM_VEC_MDPD_H
#define LMP_ATOM_VEC_MDPD_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecMDPD : public AtomVec {
 public:
  AtomVecMDPD(class LAMMPS *);
  void init();
  void force_clear(int, size_t);
  void data_atom_post(int);
  int property_atom(char *);
  void pack_property_atom(int, double *, int, int);
};

}

#endif
#endif
