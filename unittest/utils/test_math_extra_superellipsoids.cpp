/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../../src/ASPHERE/math_extra_superellipsoids.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>

// TODO: consider making a fixture with several setup functions?

static constexpr double EPSILON = 1e-4;

TEST(ContactPointAndNormal, sphere)
{
  // First grain
  double xci[3] = {1.0, 5.246, 3.123};
  double ri = 2.5;
  double shapei[3] = {ri, ri, ri};
  double Ri[3][3] = {{1.0, 0.0, 0.0},
                     {0.0, 1.0, 0.0},
                     {0.0, 0.0, 1.0}};
  double blocki[2] = {2.0, 2.0};
  int flagi = 0;

  // Second grains
  double xcj[3] = {2.0, -1.562, 4.607};
  double rj = 1.25;
  double shapej[3] = {rj, rj, rj};
  double Rj[3][3] = {{1.0, 0.0, 0.0},
                     {0.0, 1.0, 0.0},
                     {0.0, 0.0, 1.0}};
  double blockj[2] = {2.0, 2.0};
  int flagj = 0;

  // Contact detection
  double X0[4] = {0.0, 0.0, 0.0, 0.0}, nij[3];
  MathExtraSuperellipsoids::determine_contact_point(xci, Ri, shapei, blocki, flagi,
                                                    xcj, Rj, shapej, blockj, flagj,
                                                    X0, nij);
  // Analytical solution
  double X0_analytical[4] = {rj * xci[0] / (ri+rj) + ri * xcj[0] / (ri+rj),
                             rj * xci[1] / (ri+rj) + ri * xcj[1] / (ri+rj),
                             rj * xci[2] / (ri+rj) + ri * xcj[2] / (ri+rj),
                             rj / ri};
  double nij_analytical[3] = {xcj[0] - xci[0], xcj[1] - xci[1], xcj[2] - xci[2]};
  double norm = std::sqrt(nij_analytical[0]*nij_analytical[0] + nij_analytical[1]*nij_analytical[1] + nij_analytical[2]*nij_analytical[2]);
  nij_analytical[0] /= norm;
  nij_analytical[1] /= norm;
  nij_analytical[2] /= norm;

  ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON);
  ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON);
  ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON);
  ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON);

  ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
  ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
  ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);
}
