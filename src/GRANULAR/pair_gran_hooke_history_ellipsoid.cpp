/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "pair_gran_hooke_history_ellipsoid.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "math_extra.h" // probably needed for some computations
#include "math_extra_superellipsoids.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr int ITERMAX_NEWTON = 100;
static constexpr double CONVERGENCE_NEWTON = 1e-6;
static constexpr int ITERMAX_LINESEARCH = 10;
static constexpr double PARAMETER_LINESEARCH = 1e-4;
static constexpr double CUTBACK_LINESEARCH = 0.5;
static constexpr int NUMSTEP_INITIAL_GUESS = 8;

extern "C" { // General Matrices
    void dgetrf_(const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info); // Factorize
    void dgetrs_(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info); // Solve (using factorzation)
}


/* ---------------------------------------------------------------------- */

PairGranHookeHistoryEllipsoid::PairGranHookeHistoryEllipsoid(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  finitecutflag = 1;
  history = 1;
  size_history = 6;  // shear[3], previous_cp[3]

  single_extra = 10;
  svector = new double[10];

  neighprev = 0;

  nmax = 0;
  mass_rigid = nullptr;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  // keep default behavior of history[i][j] = -history[j][i]

  nondefault_history_transfer = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>(
      modify->add_fix("NEIGH_HISTORY_HH_ELL_DUMMY" + std::to_string(instance_me) + " all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryEllipsoid::~PairGranHookeHistoryEllipsoid()
{
  if (copymode) return;

  delete[] svector;

  if (!fix_history)
    modify->delete_fix("NEIGH_HISTORY_HH_ELL_DUMMY" + std::to_string(instance_me));
  else
    modify->delete_fix("NEIGH_HISTORY_HH_ELL" + std::to_string(instance_me));

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    delete[] onerad_dynamic;
    delete[] onerad_frozen;
    delete[] maxrad_dynamic;
    delete[] maxrad_frozen;
  }

  memory->destroy(mass_rigid);
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, fx, fy, fz;
  double radi, radj, radsum, rsq, r, rinv, rsqinv, factor_lj;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double wr1, wr2, wr3;
  double vtr1, vtr2, vtr3, vrel;
  double mi, mj, meff, damp, ccel, tor1, tor2, tor3;
  double fn, fs, fs1, fs2, fs3;
  double shrmag, rsht;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *shear, *allshear, **firstshear, *prev_cp; // added previous contact point placeholder

  double shapex, shapey, shapez; // ellipsoid shape params
  double quat1, quat2, quat3, quat4;
  double block1, block2;

  double X0[4], shapei[3], blocki[3], shapej[3], blockj[3], Ri[3][3], Rj[3][3];

  ev_init(eflag, vflag);

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body", tmp);
    auto *mass_body = (double *) fix_rigid->extract("masstotal", tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid, nmax, "pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0)
        mass_rigid[i] = mass_body[body[i]];
      else
        mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double *special_lj = force->special_lj;
  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec_ellipsoid->bonus;
  int *ellipsoid = atom->ellipsoid;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstshear = fix_history->firstvalue;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];


    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      radj = radius[j];
      radsum = radi + radj;

      bool touching;
      if (rsq >= radsum * radsum) {
        touching = false;
      // TODO: consider implementing a bounding-box check for hierchical detection
      //       Could be useful for high aspect ratio grain.
      //       Maybe make it an option, since it could be slower for low aspect ratio grains
      } else {
        // Super-ellipsoid contact detection between atoms i and j
        MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat, Ri);
        // TODO: Not sure if j is accessible if ghost, radius is, so bonus props must have been communicated on ghost atoms I think
        MathExtra::quat_to_mat(bonus[ellipsoid[j]].quat, Rj);

        if (touch[jj] == 1) {
          // Continued contact: use grain true shape and last contact point
          MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei);
          MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej);
          MathExtra::copy3(bonus[ellipsoid[i]].block, blocki);
          MathExtra::copy3(bonus[ellipsoid[j]].block, blockj);
          // TODO: implement neigh history!
          // TODO: move contact point with rigid body motion of the pair ?
          //       not sure if enough information to do that
          MathExtra::copy3(prev_cp, X0);
          X0[3] = 0.0; // Lagrange multiplier mu^2 initially zero
          int status = determine_contact_point(x[i], Ri, shapei, blocki, x[j], Rj, shapej, blockj, X0);
          if (status == 0)
            touching = true;
          else if(status == 5)
            touching = false;
          else
            error->all(FLERR, "Ellipsoid contact detection failed with status {} ", status);
        } else {
          // New contact: Build initial guess incrementally
          MathExtra::scaleadd3(radj / radsum, x[i], radi /radsum, x[j], X0);
          for (int iter_ig = 1 ; iter_ig <= NUMSTEP_INITIAL_GUESS ; iter_ig++) {
            X0[3] = 0.0; // Lagrange multiplier mu^2 initially zero
            double frac = iter_ig / double(NUMSTEP_INITIAL_GUESS);
            double shapei[3] = {1.0, 1.0, 1.0};
            double shapej[3] = {1.0, 1.0, 1.0};
            MathExtra::scaleadd3(1.0-frac, shapei, frac, bonus[ellipsoid[i]].shape, shapei);
            MathExtra::scaleadd3(1.0-frac, shapej, frac, bonus[ellipsoid[j]].shape, shapej);
            double blocki[2] = {2.0 + frac * (bonus[ellipsoid[i]].block[0] - 2.0), 2.0 + frac * (bonus[ellipsoid[i]].block[1] - 2.0)};
            double blockj[2] = {2.0 + frac * (bonus[ellipsoid[j]].block[0] - 2.0), 2.0 + frac * (bonus[ellipsoid[j]].block[1] - 2.0)};
            int status = determine_contact_point(x[i], Ri, shapei, blocki, x[j], Rj, shapej, blockj, X0);
            if (status == 0)
              touching = true;
            else if(status == 5)
              touching = false;
            else
              error->all(FLERR, "Ellipsoid contact detection failed with status {} ", status);
          }
        }
      }

      if (!touching) {
        // unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[3 * jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;
      } else {
        // TODO: Compute the force between the 2 superquadrics
        MathExtra::copy3(X0, prev_cp);

        // TODO: Everything below must be changed

        r = sqrt(rsq);
        rinv = 1.0 / r;
        rsqinv = 1.0 / rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1 * delx + vr2 * dely + vr3 * delz;
        vn1 = delx * vnnr * rsqinv;
        vn2 = dely * vnnr * rsqinv;
        vn3 = delz * vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        wr1 = (radi * omega[i][0] + radj * omega[j][0]) * rinv;
        wr2 = (radi * omega[i][1] + radj * omega[j][1]) * rinv;
        wr3 = (radi * omega[i][2] + radj * omega[j][2]) * rinv;

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi * mj / (mi + mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        // normal forces = Hookian contact + normal velocity damping

        damp = meff * gamman * vnnr * rsqinv;
        ccel = kn * (radsum - r) * rinv - damp;
        if (limit_damping && (ccel < 0.0)) ccel = 0.0;

        // relative velocities

        vtr1 = vt1 - (delz * wr2 - dely * wr3);
        vtr2 = vt2 - (delx * wr3 - delz * wr1);
        vtr3 = vt3 - (dely * wr1 - delx * wr2);
        vrel = vtr1 * vtr1 + vtr2 * vtr2 + vtr3 * vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;
        shear = &allshear[3 * jj];

        if (shearupdate) {
          shear[0] += vtr1 * dt;
          shear[1] += vtr2 * dt;
          shear[2] += vtr3 * dt;
        }
        shrmag = sqrt(shear[0] * shear[0] + shear[1] * shear[1] + shear[2] * shear[2]);

        if (shearupdate) {

          // rotate shear displacements

          rsht = shear[0] * delx + shear[1] * dely + shear[2] * delz;
          rsht *= rsqinv;
          shear[0] -= rsht * delx;
          shear[1] -= rsht * dely;
          shear[2] -= rsht * delz;
        }

        // tangential forces = shear + tangential velocity damping

        fs1 = -(kt * shear[0] + meff * gammat * vtr1);
        fs2 = -(kt * shear[1] + meff * gammat * vtr2);
        fs3 = -(kt * shear[2] + meff * gammat * vtr3);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1 * fs1 + fs2 * fs2 + fs3 * fs3);
        fn = xmu * fabs(ccel * r);

        if (fs > fn) {
          if (shrmag != 0.0) {
            shear[0] =
                (fn / fs) * (shear[0] + meff * gammat * vtr1 / kt) - meff * gammat * vtr1 / kt;
            shear[1] =
                (fn / fs) * (shear[1] + meff * gammat * vtr2 / kt) - meff * gammat * vtr2 / kt;
            shear[2] =
                (fn / fs) * (shear[2] + meff * gammat * vtr3 / kt) - meff * gammat * vtr3 / kt;
            fs1 *= fn / fs;
            fs2 *= fn / fs;
            fs3 *= fn / fs;
          } else
            fs1 = fs2 = fs3 = 0.0;
        }

        // forces & torques

        fx = delx * ccel + fs1;
        fy = dely * ccel + fs2;
        fz = delz * ccel + fs3;
        fx *= factor_lj;
        fy *= factor_lj;
        fz *= factor_lj;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dely * fs3 - delz * fs2);
        tor2 = rinv * (delz * fs1 - delx * fs3);
        tor3 = rinv * (delx * fs2 - dely * fs1);
        tor1 *= factor_lj;
        tor2 *= factor_lj;
        tor3 *= factor_lj;
        torque[i][0] -= radi * tor1;
        torque[i][1] -= radi * tor2;
        torque[i][2] -= radi * tor3;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= radj * tor1;
          torque[j][1] -= radj * tor2;
          torque[j][2] -= radj * tor3;
        }

        if (evflag) ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  onerad_dynamic = new double[n + 1];
  onerad_frozen = new double[n + 1];
  maxrad_dynamic = new double[n + 1];
  maxrad_frozen = new double[n + 1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::settings(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR, "Illegal pair_style command");

  kn = utils::numeric(FLERR, arg[0], false, lmp);
  if (strcmp(arg[1], "NULL") == 0)
    kt = kn * 2.0 / 7.0;
  else
    kt = utils::numeric(FLERR, arg[1], false, lmp);

  gamman = utils::numeric(FLERR, arg[2], false, lmp);
  if (strcmp(arg[3], "NULL") == 0)
    gammat = 0.5 * gamman;
  else
    gammat = utils::numeric(FLERR, arg[3], false, lmp);

  xmu = utils::numeric(FLERR, arg[4], false, lmp);
  dampflag = utils::inumeric(FLERR, arg[5], false, lmp);
  if (dampflag == 0) gammat = 0.0;

  limit_damping = 0;
  if (narg == 7) {
    if (strcmp(arg[6], "limit_damping") == 0)
      limit_damping = 1;
    else
      error->all(FLERR, "Illegal pair_style command");
  }

  if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 || xmu < 0.0 || xmu > 10000.0 ||
      dampflag < 0 || dampflag > 1)
    error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::coeff(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::init_style()
{
  int i;

  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag || !atom->angmom_flag || !atom->ellipsoid_flag)
    error->all(FLERR, "Pair gran/h/ellipsoid* requires atom attributes radius, rmass, angmom and ellipdoid flag");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Pair gran/h/ellipsoid* requires ghost atoms store velocity");

  // need a granular neighbor list

  if (history)
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_HISTORY);
  else
    neighbor->add_request(this, NeighConst::REQ_SIZE);

  dt = update->dt;

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (history && (fix_history == nullptr)) {
    auto cmd = fmt::format("NEIGH_HISTORY_HH_ELL{} all NEIGH_HISTORY {}", instance_me, size_history);
    fix_history = dynamic_cast<FixNeighHistory *>(
        modify->replace_fix("NEIGH_HISTORY_HH_ELL_DUMMY" + std::to_string(instance_me), cmd, 1));
    fix_history->pair = this;
  }

  // check for FixFreeze and set freeze_group_bit

  auto fixlist = modify->get_fix_by_style("^freeze");
  if (fixlist.size() == 0)
    freeze_group_bit = 0;
  else if (fixlist.size() > 1)
    error->all(FLERR, "Only one fix freeze command at a time allowed");
  else
    freeze_group_bit = fixlist.front()->groupbit;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = nullptr;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) {
      if (fix_rigid)
        error->all(FLERR, "Only one fix rigid command at a time allowed");
      else
        fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    for (auto &ipour : pours) {
      itype = i;
      double maxrad = *((double *) ipour->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
    for (auto &idep : deps) {
      itype = i;
      double maxrad = *((double *) idep->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
  }

  // since for ellipsoids radius is the maximum of the three axes, no need to change this part

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
  }

  MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);

  // set fix which stores history info

  if (history) {
    fix_history = dynamic_cast<FixNeighHistory *>(
        modify->get_fix_by_id("NEIGH_HISTORY_HH_ELL" + std::to_string(instance_me)));
    if (!fix_history) error->all(FLERR, "Could not find pair fix neigh history ID");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranHookeHistoryEllipsoid::init_one(int i, int j)
{
  if (!allocated) allocate();

  // cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
  cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
  cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);
  return cutoff;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) fwrite(&setflag[i][j], sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::write_restart_settings(FILE *fp)
{
  fwrite(&kn, sizeof(double), 1, fp);
  fwrite(&kt, sizeof(double), 1, fp);
  fwrite(&gamman, sizeof(double), 1, fp);
  fwrite(&gammat, sizeof(double), 1, fp);
  fwrite(&xmu, sizeof(double), 1, fp);
  fwrite(&dampflag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &kn, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &kt, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &gamman, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &gammat, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xmu, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &dampflag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&kn, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kt, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamman, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gammat, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&xmu, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dampflag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranHookeHistoryEllipsoid::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                                    double /*factor_coul*/, double /*factor_lj*/, double &fforce)
{
  double radi, radj, radsum;
  double r, rinv, rsqinv, delx, dely, delz;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3, wr1, wr2, wr3;
  double mi, mj, meff, damp, ccel;
  double vtr1, vtr2, vtr3, vrel, shrmag;
  double fs1, fs2, fs3, fs, fn;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  if (rsq >= radsum * radsum) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  r = sqrt(rsq);
  rinv = 1.0 / r;
  rsqinv = 1.0 / rsq;

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  vnnr = vr1 * delx + vr2 * dely + vr3 * delz;
  vn1 = delx * vnnr * rsqinv;
  vn2 = dely * vnnr * rsqinv;
  vn3 = delz * vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  double **omega = atom->omega;
  wr1 = (radi * omega[i][0] + radj * omega[j][0]) * rinv;
  wr2 = (radi * omega[i][1] + radj * omega[j][1]) * rinv;
  wr3 = (radi * omega[i][2] + radj * omega[j][2]) * rinv;

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  double *rmass = atom->rmass;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    // NOTE: ensure mass_rigid is current for owned+ghost atoms?
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi * mj / (mi + mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff * gamman * vnnr * rsqinv;
  ccel = kn * (radsum - r) * rinv - damp;
  if (limit_damping && (ccel < 0.0)) ccel = 0.0;

  // relative velocities

  vtr1 = vt1 - (delz * wr2 - dely * wr3);
  vtr2 = vt2 - (delx * wr3 - delz * wr1);
  vtr3 = vt3 - (dely * wr1 - delx * wr2);
  vrel = vtr1 * vtr1 + vtr2 * vtr2 + vtr3 * vtr3;
  vrel = sqrt(vrel);

  // shear history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary

  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  double *allshear = fix_history->firstvalue[i];

  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }

  double *shear = &allshear[3 * neighprev];
  shrmag = sqrt(shear[0] * shear[0] + shear[1] * shear[1] + shear[2] * shear[2]);

  // tangential forces = shear + tangential velocity damping

  fs1 = -(kt * shear[0] + meff * gammat * vtr1);
  fs2 = -(kt * shear[1] + meff * gammat * vtr2);
  fs3 = -(kt * shear[2] + meff * gammat * vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1 * fs1 + fs2 * fs2 + fs3 * fs3);
  fn = xmu * fabs(ccel * r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      fs1 *= fn / fs;
      fs2 *= fn / fs;
      fs3 *= fn / fs;
      fs *= fn / fs;
    } else
      fs1 = fs2 = fs3 = fs = 0.0;
  }

  // set force and return no energy

  fforce = ccel;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranHookeHistoryEllipsoid::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                            int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistoryEllipsoid::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairGranHookeHistoryEllipsoid::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   self-interaction range of particle
------------------------------------------------------------------------- */

double PairGranHookeHistoryEllipsoid::atom2cut(int i)
{
  double cut = atom->radius[i] * 2;
  return cut;
}

/* ----------------------------------------------------------------------
   maximum interaction range for two finite particles
------------------------------------------------------------------------- */

double PairGranHookeHistoryEllipsoid::radii2cut(double r1, double r2)
{
  double cut = r1 + r2;
  return cut;
}


// High performance versions
// TODO: this creates a fair bit of code duplication
//       but avoids recomputing some of the expensive pow(), etc that would come with creating 3 functions:
//       compute_shape(), compute_gradient, compute_jacobian.
//       not sure how to best do this without creating many small help functions
//       Pushing that logic, the calculation of a_inv, etc is not necessary. could define shapeinv
void PairGranHookeHistoryEllipsoid::derivatives_local(const double* xlocal, const double* shape, const double* block, double* grad, double hess[3][3]) {
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = xlocal[0] * a_inv;
  double y_b = xlocal[1] * b_inv;
  double z_c = xlocal[2] * c_inv;
  double n1 = block[0];
  double n2 = block[1];
  // Consider simplifying with flag_super
  double x_a_pow_n2_m2 = std::pow(std::abs(x_a), n2 - 2.0);
  double x_a_pow_n2_m1 = x_a_pow_n2_m2 * x_a;
  double y_b_pow_n2_m2 = std::pow(std::abs(y_b), n2 - 2.0);
  double y_b_pow_n2_m1 = y_b_pow_n2_m2 * y_b;

  double nu = (x_a_pow_n2_m1 * x_a) + (y_b_pow_n2_m1 * y_b);
  double nu_pow_n1_n2_m2 = std::pow(nu, n1/n2 - 2.0);
  double nu_pow_n1_n2_m1 = nu_pow_n1_n2_m1 * nu;

  double z_c_pow_n1_m2 = std::pow(std::abs(z_c), n1 -2.0);

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n1 * a_inv * x_a_pow_n2_m1 * nu_pow_n1_n2_m1 * signx;
  grad[1] = n1 * b_inv * y_b_pow_n2_m1 * nu_pow_n1_n2_m1 * signy;
  grad[2] = n1 * c_inv * (z_c_pow_n1_m2 * z_c) * signz;

  // Equation (15)
  double signxy = signx * signy;
  hess[0][0] = a_inv * a_inv * (n1 * (n2 - 1.0) * x_a_pow_n2_m2 * nu_pow_n1_n2_m1 +
                                (n1 - n2) * n1 * (x_a_pow_n2_m1 * x_a_pow_n2_m1) * nu_pow_n1_n2_m2);
  hess[1][1] = b_inv * b_inv * (n1 * (n2 - 1.0) * y_b_pow_n2_m2 * nu_pow_n1_n2_m1 +
                                (n1 - n2) * n1 * (y_b_pow_n2_m1 * y_b_pow_n2_m1) * nu_pow_n1_n2_m2);
  hess[0][1] = hess[1][0] = a_inv * b_inv * (n1 - n2) * n1 * x_a_pow_n2_m1 * y_b_pow_n2_m1 * nu_pow_n1_n2_m2 * signxy;
  hess[2][2] = c_inv * c_inv * n1 * (n1 - 1.0) * z_c_pow_n1_m2;
  hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;
}

// High performance version
double PairGranHookeHistoryEllipsoid::shape_and_gradient_local(const double* xlocal, const double* shape, const double* block, double* grad) {
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = xlocal[0] * a_inv;
  double y_b = xlocal[1] * b_inv;
  double z_c = xlocal[2] * c_inv;
  double n1 = block[0];
  double n2 = block[1];
  // Consider simplifying with flag_super
  double x_a_pow_n2_m1 = std::pow(std::abs(x_a), n2 - 1.0);
  double y_b_pow_n2_m1 = std::pow(std::abs(y_b), n2 - 1.0);

  double nu = (x_a_pow_n2_m1 * x_a) + (y_b_pow_n2_m1 * y_b);
  double nu_pow_n1_n2_m1 = std::pow(nu, n1/n2 - 1.0);

  double z_c_pow_n1_m1 = std::pow(std::abs(z_c), n1 - 1.0);

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n1 * a_inv * x_a_pow_n2_m1 * nu_pow_n1_n2_m1 * signx;
  grad[1] = n1 * b_inv * y_b_pow_n2_m1 * nu_pow_n1_n2_m1 * signy;
  grad[2] = n1 * c_inv * z_c_pow_n1_m1 * signz;

  return (nu_pow_n1_n2_m1 * nu) + (z_c_pow_n1_m1 * z_c) - 1.0;
}

double PairGranHookeHistoryEllipsoid::compute_residual(double* xci, double Ri[3][3], double* shapei, double* blocki,
                                                       double* xcj, double Rj[3][3], double* shapej, double* blockj,
                                                       double* X, double* shapefunc, double* residual) {
  double tmp[3];
  double xi_local[3], xj_local[3];
  double gradi[3], gradj[3];

  MathExtra::sub3(X, xci, tmp);
  MathExtra::transpose_matvec(Ri, tmp, xi_local);
  shapefunc[0] = shape_and_gradient_local(xi_local, shapei, blocki, tmp);
  MathExtra::matvec(Ri, tmp, gradi);

  MathExtra::sub3(X, xcj, tmp);
  MathExtra::transpose_matvec(Rj, tmp, xj_local);
  shapefunc[1] = shape_and_gradient_local(xj_local, shapej, blockj, tmp);
  MathExtra::matvec(Rj, tmp, gradj);

  // Equation (23)
  MathExtra::scaleadd3(X[3], gradj, gradi, residual);
  residual[3] = shapefunc[0] - shapefunc[1];
  return residual[0]*residual[0] + residual[1]*residual[1] + residual[2]*residual[2] + residual[3]*residual[3];
}

void PairGranHookeHistoryEllipsoid::compute_jacobian(double* xci, double Ri[3][3], double* shapei, double* blocki,
                                                     double* xcj, double Rj[3][3], double* shapej, double* blockj,
                                                     double* X, double* jacobian) {
  double tmp_v[3], tmp_m[3][3];
  double xi_local[3], xj_local[3];
  double gradi[3], hessi[3][3], gradj[3], hessj[3][3];

  MathExtra::sub3(X, xci, tmp_v);
  MathExtra::transpose_matvec(Ri, tmp_v, xi_local);
  derivatives_local(xi_local, shapei, blocki, tmp_v, hessi);
  MathExtra::matvec(Ri, tmp_v, gradi);
  MathExtra::times3_transpose(hessi, Ri, tmp_m);
  MathExtra::times3(Ri, tmp_m, hessi);

  MathExtra::sub3(X, xcj, tmp_v);
  MathExtra::transpose_matvec(Rj, tmp_v, xj_local);
  derivatives_local(xj_local, shapej, blockj, tmp_v, hessj);
  MathExtra::matvec(Rj, tmp_v, gradj);
  MathExtra::times3_transpose(hessj, Rj, tmp_m);
  MathExtra::times3(Rj, tmp_m, hessj);

  // Jacobian (derivative of residual)
  // 1D column-major matrix for LAPACK/linalg compatibility
  for (int row = 0 ; row < 3 ; row++) {
    for (int col = 0 ; col < 3 ; col++) {
      jacobian[row + col*4] = hessi[row][col] + X[3] * hessj[row][col];
    }
    jacobian[row + 3*4] = gradj[row];
  }
  for (int col = 0 ; col < 3 ; col++) {
    jacobian[3 + col*4] = gradi[col] - gradj[col];
  }
  jacobian[15] = 0.0;
}


int PairGranHookeHistoryEllipsoid::determine_contact_point(double* xci, double Ri[3][3], double* shapei, double* blocki,
                                                            double* xcj, double Rj[3][3], double* shapej, double* blockj,
                                                            double* X0) {
  double norm, norm_ini, shapefunc[2], residual[4], jacobian[16];
  bool converged(false);
  norm = compute_residual(xci, Ri, shapei, blocki, xcj, Rj, shapej, blockj, X0, shapefunc, residual);
  for (int iter = 0 ; iter < ITERMAX_NEWTON ; iter++) {
    norm_ini = norm;
    compute_jacobian(xci, Ri, shapei, blocki, xcj, Rj, shapej, blockj, X0, jacobian);

    // Solve Newton step
    int lapack_error, ipiv[16];
    const int n = 4;
    const char trans = 'N';
    const int nrhs = 1;
    double rhs[4] = {-residual[0], -residual[1], -residual[2], -residual[3]};
    dgetrf_(&n, &n, jacobian, &n, ipiv, &lapack_error);
    if (lapack_error)
      return lapack_error;
    dgetrs_(&trans, &n, &nrhs, jacobian, &n, ipiv, rhs, &n, &lapack_error);
    if (lapack_error)
      return lapack_error;

    // Backtracking line search
    double a(1.0), X_line[4];
    for (int iter_ls = 0 ; iter_ls < ITERMAX_LINESEARCH ; iter_ls++) {
      X_line[0] = X0[0] + a * rhs[0];
      X_line[1] = X0[1] + a * rhs[1];
      X_line[2] = X0[2] + a * rhs[2];
      X_line[3] = X0[3] + a * rhs[3];

      norm = compute_residual(xci, Ri, shapei, blocki, xcj, Rj, shapej, blockj, X_line, shapefunc, residual);
      if (norm < norm_ini - PARAMETER_LINESEARCH * a * norm_ini)
        break; // Armijo - Goldstein condition
      else
        a *= CUTBACK_LINESEARCH;
    }
    X0[0] = X_line[0];
    X0[1] = X_line[1];
    X0[2] = X_line[2];
    X0[3] = X_line[3];

    if (norm < CONVERGENCE_NEWTON) {
      converged = true;
      break;
    }
  }

  // LAPACK error are within [-4, 4], use 5 non-touching, -5 non-converging
  if (!converged)
    return -5;
  if (shapefunc[0] > 0.0 || shapefunc[1] > 0.0)
    return 5;

  return 0;
}
