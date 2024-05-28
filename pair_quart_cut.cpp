
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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_quart_cut.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairQUARTCut::PairQUARTCut(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairQUARTCut::~PairQUARTCut()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(prefactor0);
    memory->destroy(prefactor1);
    memory->destroy(prefactor2);
    memory->destroy(quart0);
    memory->destroy(quart1);
    memory->destroy(quart2);
    memory->destroy(quart3);
    memory->destroy(quart4);
    memory->destroy(quart5);
    memory->destroy(quart6);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairQUARTCut::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, drr0, forcequart, factor_quart;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_quart = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        drr0 = sqrt(rsq)-quart0[itype][jtype];
        forcequart = - 1.0 * drr0 * (quart1[itype][jtype] + quart2[itype][jtype] * drr0 + quart3[itype][jtype] * pow(drr0,2))/sqrt(rsq);
        fpair = factor_quart * forcequart;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          evdwl = drr0 * drr0 * (quart4[itype][jtype] + quart5[itype][jtype]*drr0 + quart6[itype][jtype] * pow(drr0,2.0)) + prefactor0[itype][jtype] - offset[itype][jtype];
          evdwl *= factor_quart;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQUARTCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(prefactor0, n, n, "pair:prefactor0");
  memory->create(prefactor1, n, n, "pair:prefactor1");
  memory->create(prefactor2, n, n, "pair:prefactor2");
  memory->create(quart0, n, n, "pair:quart0");
  memory->create(quart1, n, n, "pair:quart1");
  memory->create(quart3, n, n, "pair:quart3");
  memory->create(quart4, n, n, "pair:quart4");
  memory->create(quart5, n, n, "pair:quart5");
  memory->create(quart6, n, n, "pair:quart6");
  memory->create(offset, n, n, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQUARTCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQUARTCut::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double prefactor0_one = utils::numeric(FLERR, arg[4], false, lmp);
  double prefactor1_one = utils::numeric(FLERR, arg[5], false, lmp);
  double prefactor2_one = utils::numeric(FLERR, arg[6], false, lmp);
  double cut_one = cut_global;
  if (narg == 8) cut_one = utils::numeric(FLERR, arg[7], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      prefactor0[i][j] = prefactor0_one;
      prefactor1[i][j] = prefactor1_one;
      prefactor2[i][j] = prefactor2_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairQUARTCut::init_style()
{
  // request regular or rRESPA neighbor list

  int list_style = NeighConst::REQ_DEFAULT;
  
  neighbor->add_request(this, list_style);

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQUARTCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    prefactor0[i][j] = 0.5 * (prefactor0[i][j] + prefactor0[j][i]);
    prefactor1[i][j] = 0.5 * (prefactor1[i][j] + prefactor1[j][i]);
    prefactor2[i][j] = 0.5 * (prefactor2[i][j] + prefactor2[j][i]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  quart0[i][j] = pow(2.0,1.0/6.0) * sigma[i][j];
  quart1[i][j] = 36.0 * pow(4.0, 1.0/3.0) * epsilon[i][j] * pow(sigma[i][j],-2.0);
  quart2[i][j] = -378.0 * prefactor1[i][j] * sqrt(2.0) * epsilon[i][j] * pow(sigma[i][j],-3.0);
  quart3[i][j] = 2226.0 *  prefactor2[i][j]* pow(2.0, 1.0/3.0) * epsilon[i][j] * pow(sigma[i][j],-4.0);

  quart4[i][j] = 0.5 * 36.0  * pow(4.0, 1.0/3.0) * epsilon[i][j] * pow(sigma[i][j],-2.0);
  quart5[i][j] = -378.0 * 1.0/3.0 * prefactor1[i][j] * sqrt(2.0) * epsilon[i][j] * pow(sigma[i][j],-3.0);
  quart6[i][j] = 2226.0 * 0.25 * prefactor2[i][j] * pow(2.0, 1.0/3.0) * epsilon[i][j] * pow(sigma[i][j],-4.0);

 

  if (offset_flag && (cut[i][j] > 0.0)) {
    double diff = cut[i][j] - quart0[i][j]  ;
    offset[i][j] = quart4[i][j] * pow(diff,2.0) + quart5[i][j] * pow(diff,3.0) + quart6[i][j] * pow(diff,4.0) + prefactor0[i][j];
  } else
    offset[i][j] = 0.0;

  quart0[j][i] = quart0[i][j];
  quart1[j][i] = quart1[i][j];
  quart2[j][i] = quart2[i][j];
  quart3[j][i] = quart3[i][j];
  quart4[j][i] = quart4[i][j];
  quart5[j][i] = quart5[i][j];
  quart6[j][i] = quart6[i][j];
  prefactor0[j][i] =  prefactor0[i][j];
  prefactor1[j][i] =  prefactor1[i][j];
  prefactor2[j][i] =  prefactor2[i][j];
  offset[j][i] = offset[i][j];



  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQUARTCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&prefactor0[i][j], sizeof(double), 1, fp);
        fwrite(&prefactor1[i][j], sizeof(double), 1, fp);
        fwrite(&prefactor2[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQUARTCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &prefactor0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &prefactor1[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &prefactor2[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&prefactor0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&prefactor1[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&prefactor2[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQUARTCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQUARTCut::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairQUARTCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g %g %g %g\n", i, epsilon[i][i], sigma[i][i], prefactor0[i][i], prefactor1[i][i], prefactor2[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairQUARTCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], prefactor0[i][j], prefactor1[i][j], prefactor2[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */
double PairQUARTCut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                            double /*factor_coul*/, double factor_quart,
                            double &fforce)
{
  double drr0, philj, forcequart;

  drr0 = sqrt(rsq)-quart0[itype][jtype];
  forcequart = - 1.0 * drr0 * (quart1[itype][jtype] + quart2[itype][jtype] * drr0 + quart3[itype][jtype] * pow(drr0,2))/sqrt(rsq);
  fforce = factor_quart * forcequart;
  philj = drr0 * drr0 * (quart4[itype][jtype] + quart5[itype][jtype]*drr0 + quart6[itype][jtype] * pow(drr0,2.0)) + prefactor0[itype][jtype] - offset[itype][jtype];
  return factor_quart * philj;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void *PairQUARTCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  if (strcmp(str, "prefactor0") == 0) return (void *) prefactor0;
  if (strcmp(str, "prefactor1") == 0) return (void *) prefactor1;
  if (strcmp(str, "prefactor2") == 0) return (void *) prefactor2;
  return nullptr;
}
