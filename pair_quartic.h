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

#ifdef PAIR_CLASS

PairStyle(quartic,PairQUARTIC)

#else

#ifndef LMP_PAIR_QUARTIC_H
#define LMP_PAIR_QUARTIC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairQUARTIC : public Pair {
 public:
  PairQUARTIC(class LAMMPS *);
  virtual ~PairQUARTIC();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  void write_data(FILE *QUARTIC);
  void write_data_all(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global;
  double **cut;
  double **quar_A0,**quar_A2,**quar_A3,**quar_A4,**quar_r0;
  double **lj1,**lj2,**lj3,**lj4;
  double **ljsw1,**ljsw2,**ljsw3,**ljsw4,**ljsw5;
  

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/