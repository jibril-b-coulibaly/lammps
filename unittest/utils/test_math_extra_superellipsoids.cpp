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
#include "../../src/math_extra.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>
#include <vector>

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
  MathExtra::norm3(nij_analytical);

  ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON);
  ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON);
  ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON);
  ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON);

  ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
  ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
  ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);

  // Rotational invariance
  double anglei = 0.456;
  double axisi[3] = {1,2,3};
  MathExtra::norm3(axisi);
  double quati[4] = {std::cos(anglei),
                     std::sin(anglei)*axisi[0],
                     std::sin(anglei)*axisi[1],
                     std::sin(anglei)*axisi[2]};
  MathExtra::quat_to_mat(quati, Ri);

  double anglej = 0.123;
  double axisj[3] = {-1,2,1};
  MathExtra::norm3(axisj);
  double quatj[4] = {std::cos(anglej),
                     std::sin(anglej)*axisj[0],
                     std::sin(anglej)*axisj[1],
                     std::sin(anglej)*axisj[2]};
  MathExtra::quat_to_mat(quatj, Rj);

  X0[0] = X0[1] = X0[2] = X0[3] = 0.0;
  MathExtraSuperellipsoids::determine_contact_point(xci, Ri, shapei, blocki, flagi,
                                                    xcj, Rj, shapej, blockj, flagj,
                                                    X0, nij);

  ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON);
  ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON);
  ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON);
  ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON);

  ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
  ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
  ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);

 
}

TEST(ContactPointAndNormal, supersphere_mono)
{
  double r = 3.456;
  double xci[3] = {-2*r, 0.0, 0.0};
  double xcj[3] = {2*r, 0.0, 0.0};
  double shape[3] = {r, r, r};
  double R[3][3] = {{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}};

  std::vector<double> blocks = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
  for (auto n : blocks) {
    double block[2] = {n, n};
    int flag = MathExtraSuperellipsoids::determine_flag(block);

    // Contact detection
    // Some starting point away from (0,0,0). Possibly bad initial guess so test is demanding
    double X0[4] = {r, -r, 2*r, 0.0}, nij[3];
    int status = MathExtraSuperellipsoids::determine_contact_point(xci, R, shape, block, flag,
                                                                   xcj, R, shape, block, flag,
                                                                   X0, nij);
    // Analytical solution
    double X0_analytical[4] = {0.0, 0.0, 0.0, 1.0};
    double nij_analytical[3] = {1.0, 0.0, 0.0};
    // TODO / WIP:
    // Gradients can be smaller in different directions, hard to naviguate canyon on high blockiness
    // Little progress made along the flat faces. Maybe use Levenberg-Marquardt or Newton with momentum (previous step memory) or other methods.
    // I think this might be an intrinsic problem with the solution having a Hessian of zero, leading to slow convergence + bad conditioning in Newton's method.

    std::cout<<n<<" "<<status<<" "<<X0[0]<<" "<<X0[1]<<" "<<X0[2]<<" "<<X0[3]<<std::endl;
    ASSERT_NEAR(X0[0], X0_analytical[0], EPSILON);
    ASSERT_NEAR(X0[1], X0_analytical[1], EPSILON);
    ASSERT_NEAR(X0[2], X0_analytical[2], EPSILON);
    ASSERT_NEAR(X0[3], X0_analytical[3], EPSILON);

    ASSERT_NEAR(nij[0], nij_analytical[0], EPSILON);
    ASSERT_NEAR(nij[1], nij_analytical[1], EPSILON);
    ASSERT_NEAR(nij[2], nij_analytical[2], EPSILON);
  }
}


// for polydisperse solution should be at the radii ratio