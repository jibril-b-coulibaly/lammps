/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_viscous_rigid.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "fix_rigid_small.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixViscousRigid::FixViscousRigid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  idrigid(NULL),fixrigid(NULL),gamma_linear(0),gamma_angular(0)
{
  dynamic_group_allow = 1;

  if (narg < 6) error->all(FLERR,"Illegal fix viscous/rigid command");
  
  int n = strlen(arg[3]) + 1;
  idrigid = new char[n];
  strcpy(idrigid,arg[3]);

  gamma_linear = force->numeric(FLERR,arg[4]);
  gamma_angular = force->numeric(FLERR,arg[5]);

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixViscousRigid::~FixViscousRigid()
{
  delete [] idrigid;
}

/* ---------------------------------------------------------------------- */

int FixViscousRigid::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscousRigid::init()
{
  int max_respa = 0;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = max_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,max_respa);
  }
	
	// set fixrigid
  
  fixrigid = NULL;
  int ifix = modify->find_fix(idrigid);
	if (ifix < 0)
		error->all(FLERR,"FixRigidSmall ID for fix viscous/rigid does not exist");
	fixrigid = (FixRigidSmall *) modify->fix[ifix];

	int flag = 0;
	if (strstr(fixrigid->style,"rigid/") == NULL) flag = 1;
	if (strstr(fixrigid->style,"/small") == NULL) flag = 1;
	if (flag)
    error->all(FLERR,"Fix viscous/rigid does not use fix rigid/small fix");

	for (int j = 0; j < ifix; j++)
	  if(strcmp(modify->fix[j]->style,"viscous/rigid") == 0)
	    error->warning(FLERR,"Fix viscous/rigid must be defined after fix rigid");
	
	if (!fixrigid->earlyflag) {
	  char str[128];
	  snprintf(str,128,"Fix %s alters rigid body forces before fix rigid defines them: useless",this->id);
	  error->warning(FLERR,str);
	} else {
	  char str[256];
	  snprintf(str,256,"Fix %s rightfully alters rigid body forces after fix rigid. Disregard warnings from fix rigid/small for this fix",this->id);
	  error->warning(FLERR,str);
	}
}

/* ---------------------------------------------------------------------- */

void FixViscousRigid::setup(int vflag)
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

void FixViscousRigid::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousRigid::post_force(int /*vflag*/)
{
  // apply drag force to rigid bodies, direction is opposed to velocity vector
	// apply drag torque to rigid bodies, direction is opposed to angular velocity vector
	// viscous forces and torques directly applied to rigid bodies, not to atoms

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
	
	int ibody;
	FixRigidSmall::Body *body;
	for (int i = 0; i < nlocal; i++) {
		if (!(mask[i] & groupbit)) continue;
	  ibody = fixrigid->bodyown[i];
	  if (ibody < 0) continue;
	  body = &fixrigid->body[ibody];
		
		body->fcm[0] -= gamma_linear*body->vcm[0];
		body->fcm[1] -= gamma_linear*body->vcm[1];
		body->fcm[2] -= gamma_linear*body->vcm[2];
		
		body->torque[0] -= gamma_angular*body->omega[0];
		body->torque[1] -= gamma_angular*body->omega[1];
		body->torque[2] -= gamma_angular*body->omega[2];
	}
}

/* ---------------------------------------------------------------------- */

void FixViscousRigid::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousRigid::min_post_force(int vflag)
{
  post_force(vflag);
}
