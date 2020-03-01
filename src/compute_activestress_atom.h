/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMPUTE_CLASS

ComputeStyle(activestress/atom,ComputeActivestressAtom)

#else

#ifndef LMP_COMPUTE_ACTIVESTRESS_ATOM_H
#define LMP_COMPUTE_ACTIVESTRESS_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeActivestressAtom : public Compute {
 public:
  ComputeActivestressAtom(class LAMMPS *, int, char **);
  ~ComputeActivestressAtom();
  void init() {}
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag,biasflag;
  Compute *temperature;
  char *id_temp;

  int nmax;
  double **stress, f_active;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute stress/atom temperature ID

Self-explanatory.

E: Compute stress/atom temperature ID does not compute temperature

The specified compute must compute temperature.

E: Per-atom virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to have
tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
