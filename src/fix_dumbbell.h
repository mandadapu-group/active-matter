/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of custom LAMMPS code for simulations of active matter.
    author: Cory Hargus
    e-mail: hargus@berkeley.edu
    github: https://github.com/mandadapu-group/active-matter
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(dumbbell,FixDumbbell)

#else

#ifndef LMP_FIX_DUMBBELL_H
#define LMP_FIX_DUMBBELL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDumbbell : public Fix {
 public:
  FixDumbbell(class LAMMPS *, int, char **);
  virtual ~FixDumbbell() {}
  int setmask();
  virtual void post_force(int);

 protected:
  class RanPark *random;
  int activestyle;
  double f_active;
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
