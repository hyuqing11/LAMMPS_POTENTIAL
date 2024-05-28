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

/* ----------------------------------------------------------------------
   Contributing author: Mark Stevens (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_quartic.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairQUARTIC::PairQUARTIC(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairQUARTIC::~PairQUARTIC()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    //memory->destroy(cut_inner_sq);
    memory->destroy(quar_A0);
    memory->destroy(quar_A2);
    memory->destroy(quar_A3);
    memory->destroy(quar_A4);
    memory->destroy(quar_r0);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(ljsw1);
    memory->destroy(ljsw2);
    memory->destroy(ljsw3);
    memory->destroy(ljsw4);
    memory->destroy(ljsw5);
  }
}

/* ---------------------------------------------------------------------- */

void PairQUARTIC::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,drVal;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  double r,t,fswitch,eswitch;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

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
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        //r2inv = 1.0/rsq;
        //r6inv = r2inv*r2inv*r2inv;
        //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        r = sqrt(rsq);
        drVal = r - quar_r0[itype][jtype];
        forcelj= -1.0*(2*quar_A2[itype][jtype]*drVal+3*quar_A3[itype][jtype]*pow(drVal,2)+4*quar_A4[itype][jtype]*pow(drVal,3))/r;
        
        /*if (rsq > cut_inner_sq[itype][jtype]) {
          //r = sqrt(rsq);
          t = r - cut_inner[itype][jtype];
          fswitch = r*t*t*(ljsw1[itype][jtype] + ljsw2[itype][jtype]*t);
          forcelj += fswitch;
        }*/
        fpair = forcelj;

       // fpair = factor_lj*forcelj*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = quar_A0[itype][jtype]+quar_A2[itype][jtype]*pow(drVal,2)+quar_A3[itype][jtype]*pow(drVal,3)+quar_A4[itype][jtype]*pow(drVal,4);
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQUARTIC::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  //memory->create(cut_inner,n+1,n+1,"pair:cut_inner");
  //memory->create(cut_inner_sq,n+1,n+1,"pair:cut_inner_sq");
  memory->create(quar_A0,n+1,n+1,"pair:quar_A0");
  memory->create(quar_A2,n+1,n+1,"pair:quar_A2");
  memory->create(quar_A3,n+1,n+1,"pair:quar_A3");
  memory->create(quar_A4,n+1,n+1,"pair:quar_A4");
   memory->create(quar_r0,n+1,n+1,"pair:quar_r0");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(ljsw1,n+1,n+1,"pair:ljsw1");
  memory->create(ljsw2,n+1,n+1,"pair:ljsw2");
  memory->create(ljsw3,n+1,n+1,"pair:ljsw3");
  memory->create(ljsw4,n+1,n+1,"pair:ljsw4");
  memory->create(ljsw5,n+1,n+1,"pair:ljsw5");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQUARTIC::settings(int narg, char **arg)
{
  if (narg != 1) {
   
    error->all(FLERR,"Illegal pair_style command");
}
   

  //cut_inner_global = force->numeric(FLERR,arg[0]);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);

  //if (cut_inner_global > cut_global)
   // error->all(FLERR,"Illegal pair_style command");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          //cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQUARTIC::coeff(int narg, char **arg)
{
  if (narg != 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double quar_A0_one= utils::numeric(FLERR,arg[2],false,lmp);
  double quar_A2_one = utils::numeric(FLERR,arg[3],false,lmp);
  double quar_A3_one = utils::numeric(FLERR,arg[4],false,lmp);
  double quar_A4_one = utils::numeric(FLERR,arg[5],false,lmp);
  double quar_r0_one = utils::numeric(FLERR,arg[6],false,lmp);
  double cut_one = utils::numeric(FLERR,arg[7],false,lmp);
  //double cut_inner_one = cut_inner_global;
 /* double cut_one = cut_global;
  if (narg == 6) {
    cut_inner_one = utils::numeric(FLERR,arg[4],false,lmp);
    cut_one = utils::numeric(FLERR,arg[5],false,lmp);
  }

  if (cut_inner_one <= 0.0 || cut_inner_one > cut_one)
    error->all(FLERR,"Incorrect args for pair coefficients");*/

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      quar_A0[i][j] = quar_A0_one;
      quar_A2[i][j] = quar_A2_one;
      quar_A3[i][j] = quar_A3_one;
      quar_A4[i][j] = quar_A4_one;
      quar_r0[i][j] = quar_r0_one;
      //cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQUARTIC::init_one(int i, int j)
{
  /*if (setflag[i][j] == 0) {
    quar_A0[i][j] = mix_energy(quar_A0[i][i],quar_A0[j][j],
                               quar_A0[i][i],quar_A0[j][j]);
    quar_A2[i][j] = mix_distance(quar_A2[i][i],quar_A2[j][j]);
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }*/

  //cut_inner_sq[i][j] = cut_inner[i][j]*cut_inner[i][j];
  quar_A0[j][i] = quar_A0[i][j];
  quar_A2[j][i] = quar_A2[i][j];
  quar_A3[j][i] = quar_A3[i][j];
  quar_A4[j][i] = quar_A4[i][j];
  quar_r0[j][i] = quar_r0[i][j];
  /*lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  double r6inv = 1.0/pow(cut[i][j],6.0);
  double r8inv = 1.0/pow(cut[i][j],8.0);
  double t = cut[i][j] - cut_inner[i][j];
  double t2inv = 1.0/(t*t);
  double t3inv = t2inv/t;
  double t3 = 1.0/t3inv;
  double a6 = (7.0*cut_inner[i][j] - 10.0*cut[i][j])*r8inv*t2inv;
  double b6 = (9.0*cut[i][j] -  7.0*cut_inner[i][j])*r8inv*t3inv;
  double a12 = (13.0*cut_inner[i][j] - 16.0*cut[i][j])*r6inv*r8inv*t2inv;
  double b12 = (15.0*cut[i][j] - 13.0*cut_inner[i][j])*r6inv*r8inv*t3inv;
  double c6 = r6inv - t3*(6.0*a6/3.0 + 6.0*b6*t/4.0);
  double c12 = r6inv*r6inv - t3*(12.0*a12/3.0 + 12.0*b12*t/4.0);

  ljsw1[i][j] = lj1[i][j]*a12 - lj2[i][j]*a6;
  ljsw2[i][j] = lj1[i][j]*b12 - lj2[i][j]*b6;
  ljsw3[i][j] = -lj3[i][j]*12.0*a12/3.0 + lj4[i][j]*6.0*a6/3.0;
  ljsw4[i][j] = -lj3[i][j]*12.0*b12/4.0 + lj4[i][j]*6.0*b6/4.0;
  ljsw5[i][j] = -lj3[i][j]*c12 + lj4[i][j]*c6;

  //cut_inner[j][i] = cut_inner[i][j];
  //cut_inner_sq[j][i] = cut_inner_sq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  ljsw1[j][i] = ljsw1[i][j];
  ljsw2[j][i] = ljsw2[i][j];
  ljsw3[j][i] = ljsw3[i][j];
  ljsw4[j][i] = ljsw4[i][j];
  ljsw5[j][i] = ljsw5[i][j];*/

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQUARTIC::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&quar_A0[i][j],sizeof(double),1,fp);
        fwrite(&quar_A2[i][j],sizeof(double),1,fp);
        fwrite(&quar_A3[i][j],sizeof(double),1,fp);
        fwrite(&quar_A4[i][j],sizeof(double),1,fp);
        fwrite(&quar_r0[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQUARTIC::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&quar_A0[i][j],sizeof(double),1,fp);
          fread(&quar_A2[i][j],sizeof(double),1,fp);
          fread(&quar_A3[i][j],sizeof(double),1,fp);
          fread(&quar_A4[i][j],sizeof(double),1,fp);
          fread(&quar_r0[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&quar_A0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&quar_A2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&quar_A3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&quar_A4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&quar_r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQUARTIC::write_restart_settings(FILE *fp)
{
  //fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQUARTIC::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    //fread(&cut_inner_global,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  //MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairQUARTIC::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,quar_A0[i][i],quar_A2[i][i],quar_A3[i][i],quar_A4[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairQUARTIC::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
              quar_A0[i][i],quar_A2[i][i],quar_A3[i][i],quar_A4[i][i],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairQUARTIC::single(int /*i*/, int /*j*/, int itype, int jtype,
                             double rsq,
                             double /*factor_coul*/, double factor_lj,
                             double &fforce)
{
  double r2inv,r6inv,forcelj,philj,drVal;
  double r,t,fswitch,phiswitch;

  //r2inv = 1.0/rsq;
  //r6inv = r2inv*r2inv*r2inv;
  r = sqrt(rsq);
  drVal = r - quar_r0[itype][jtype];
  forcelj= -1.0*(2*quar_A2[itype][jtype]*drVal+3*quar_A3[itype][jtype]*pow(drVal,2)+4*quar_A4[itype][jtype]*pow(drVal,3))/r;
  //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  /*if (rsq > cut_inner_sq[itype][jtype]) {
    r = sqrt(rsq);
    t = r - cut_inner[itype][jtype];
    fswitch = r*t*t*(ljsw1[itype][jtype] + ljsw2[itype][jtype]*t);
    forcelj += fswitch;
  }*/
  //fforce = factor_lj*forcelj*r2inv;
  fforce = forcelj; 
  //philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
  //philj += ljsw5[itype][jtype];
  /*if (rsq > cut_inner_sq[itype][jtype]) {
    phiswitch = t*t*t*(ljsw3[itype][jtype] + ljsw4[itype][jtype]*t);
    philj += phiswitch;
  }*/

  return fforce;
}