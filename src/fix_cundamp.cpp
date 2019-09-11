/* ----------------------------------------------------------------------
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

#include "fix_cundamp.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCundamp::FixCundamp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  gamma_linear(NULL),gamma_angular(NULL)
{
  dynamic_group_allow = 1;
  
  if (!atom->sphere_flag)
	error->all(FLERR,"Fix cundamp requires atom style sphere");

  if (narg < 5) error->all(FLERR,"Illegal fix cundamp command");
	
	double gamma_linear_one = force->numeric(FLERR,arg[3]);
	double gamma_angular_one = force->numeric(FLERR,arg[4]);
  gamma_linear = new double[atom->ntypes+1];
	gamma_angular = new double[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) {
		gamma_linear[i] = gamma_linear_one;
		gamma_angular[i] = gamma_angular_one;
	}

  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix cundamp command");
      int itype = force->inumeric(FLERR,arg[iarg+1]);
      double scale = force->numeric(FLERR,arg[iarg+2]);
      if (itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Illegal fix cundamp command");
      gamma_linear[itype] = gamma_linear_one * scale;
			gamma_angular[itype] = gamma_angular_one * scale;
      iarg += 3;
    } else error->all(FLERR,"Illegal fix cundamp command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixCundamp::~FixCundamp()
{
  delete [] gamma_linear;
	delete [] gamma_angular;
}

/* ---------------------------------------------------------------------- */

int FixCundamp::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCundamp::init()
{
  int max_respa = 0;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = max_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,max_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixCundamp::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixCundamp::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCundamp::post_force(int /*vflag*/)
{
  // apply force reduction/increase to granular particles. Force is reduced/increased if its power is positive/negative
  // apply torque reduction/increase to granular partiicles. Torque is reduced/increased if its power is positive/negative
  // applied over each direction independently -> artificial, non-objective, frame-dependent damping method
  
  
  double **v = atom->v;
	double **omega = atom->omega;
	double **f = atom->f;
	double **torque = atom->torque;
	int *mask = atom->mask;
	int *type = atom->type;
  int nlocal = atom->nlocal;
	
	double gamma_l,gamma_a;
	int sgnf0,sgnf1,sgnf2;
	int sgntq0,sgntq1,sgntq2;
	
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
			gamma_l = gamma_linear[type[i]];
			gamma_a = gamma_angular[type[i]];
			
			sgnf0 = (f[i][0]*v[i][0] > 0) - (f[i][0]*v[i][0] < 0);
			sgnf1 = (f[i][1]*v[i][1] > 0) - (f[i][1]*v[i][1] < 0);
			sgnf2 = (f[i][2]*v[i][2] > 0) - (f[i][2]*v[i][2] < 0);
			f[i][0] *= 1.0-gamma_l*sgnf0;
			f[i][1] *= 1.0-gamma_l*sgnf1;
			f[i][2] *= 1.0-gamma_l*sgnf2;
			
			sgntq0 = (torque[i][0]*omega[i][0] > 0) - (torque[i][0]*omega[i][0] < 0);
			sgntq1 = (torque[i][1]*omega[i][1] > 0) - (torque[i][1]*omega[i][1] < 0);
			sgntq2 = (torque[i][2]*omega[i][2] > 0) - (torque[i][2]*omega[i][2] < 0);
			torque[i][0] *= 1.0-gamma_a*sgntq0;
			torque[i][1] *= 1.0-gamma_a*sgntq1;
			torque[i][2] *= 1.0-gamma_a*sgntq2;
    }
}

/* ---------------------------------------------------------------------- */

void FixCundamp::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCundamp::min_post_force(int vflag)
{
  post_force(vflag);
}
