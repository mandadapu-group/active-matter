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

#include <mpi.h>
#include <cstring>
#include "compute_activestress.h"
#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include <iostream>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------

Compute the 2x2 active stress tensor for a 2D system of torqued dumbbells.
Assumes counter-clockwise active torque.

---------------------------------------------------------------------- */

ComputeActiveStress::ComputeActiveStress(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal compute activestress command");
  f_active = force->numeric(FLERR,arg[3]);
  int dimension = domain->dimension;
  if (dimension != 2) error->all(FLERR,"Only 2D stress tensor is supported.");

  vector_flag = 1;
  size_vector = 4;      // unraveled 2x2 active stress tensor
  vector = new double[4];
}

/* ---------------------------------------------------------------------- */

ComputeActiveStress::~ComputeActiveStress()
{
  delete [] vector;
}

void ComputeActiveStress::active_compute()
{
  int i1, i2;
  double T_A[4];
  double delx, dely, rsq, r;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  for (int n = 0; n < 4; n++) T_A[n] = 0.0;
  for (int n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    factor = 0.0;
    if (i1 < nlocal)
      factor += 0.5;
    if (i2 < nlocal)
      factor += 0.5;
    delx = x[i2][0] - x[i1][0];
    dely = x[i2][1] - x[i1][1];
    rsq = delx*delx + dely*dely;
    r = sqrt(rsq);
    T_A[0] += factor * f_active * dely * delx / r;
    T_A[1] += factor * f_active * dely * dely / r;
    T_A[2] -= factor * f_active * delx * delx / r;
    T_A[3] -= factor * f_active * delx * dely / r;
  }
  MPI_Allreduce(T_A,vector,4,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void ComputeActiveStress::compute_vector()
{
  invoked_vector = update->ntimestep;
  nktv2p = force->nktv2p;
  inv_volume = 1.0 / (domain->xprd * domain->yprd);
  active_compute();
  for (int i = 0; i < 4; i++)
    vector[i] *= inv_volume * nktv2p;
}
