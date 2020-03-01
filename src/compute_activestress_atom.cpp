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

#include <cstdlib>
#include <cstring>
#include <cmath>
#include "compute_activestress_atom.h"
#include "neighbor.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------

Compute the per-atom contributions to the 2x2 active stress tensor for a
2D system of torqued dumbbells. Assumes counter-clockwise active torque.

---------------------------------------------------------------------- */

ComputeActivestressAtom::ComputeActivestressAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_temp(NULL), stress(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute activestress/atom command");
  f_active = force->numeric(FLERR,arg[3]);

  peratom_flag = 1;
  size_peratom_cols = 4;
  pressatomflag = 1;
  timeflag = 1;
  comm_reverse = 4;
  nmax = 0;
  keflag = 0;
  pairflag = 0;
  bondflag = 0;
}

/* ---------------------------------------------------------------------- */

ComputeActivestressAtom::~ComputeActivestressAtom()
{
  delete [] id_temp;
  memory->destroy(stress);
}

/* ---------------------------------------------------------------------- */

void ComputeActivestressAtom::compute_peratom()
{
  int i,j;
  double onemass;

  invoked_peratom = update->ntimestep;
  if (update->vflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom virial was not tallied on needed timestep");

  // grow local stress array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(stress);
    nmax = atom->nmax;
    memory->create(stress,nmax,4,"stress/atom:stress");
    array_atom = stress;
  }

  // ntotal includes ghosts if either newton flag is set
  int nlocal = atom->nlocal;
  int ntotal = nlocal;
  if (force->newton) ntotal += atom->nghost;

  // clear local stress array
  for (i = 0; i < nlocal; i++)
    for (j = 0; j < 4; j++)
      stress[i][j] = 0.0;

  int *mask = atom->mask;
  // include the asymmetric active force contribution to the stress tensor
  double delx, dely, rsq, r;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **x = atom->x;
  for (i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    // Get a unit vector pointing from atom 1 to atom 2 (assuming 2d in xy-plane)
    delx = x[i2][0] - x[i1][0];
    dely = x[i2][1] - x[i1][1];
    rsq = delx*delx + dely*dely;
    r = sqrt(rsq);

    if (i1 < nlocal)
      stress[i1][0] += .5 * f_active * dely * delx / r;
      stress[i1][1] += .5 * f_active * dely * dely / r;
      stress[i1][2] += -.5 * f_active * delx * delx / r;
      stress[i1][3] += -.5 * f_active * delx * dely / r;

    if (i2 < nlocal)
      stress[i2][0] += .5 * f_active * dely * delx / r;
      stress[i2][1] += .5 * f_active * dely * dely / r;
      stress[i2][2] += -.5 * f_active * delx * delx / r;
      stress[i2][3] += -.5 * f_active * delx * dely / r;
  }

  // convert to stress*volume units = -pressure*volume
  double nktv2p = -force->nktv2p;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      stress[i][0] *= nktv2p;
      stress[i][1] *= nktv2p;
      stress[i][2] *= nktv2p;
      stress[i][3] *= nktv2p;
    }
}

/* ---------------------------------------------------------------------- */

int ComputeActivestressAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = stress[i][0];
    buf[m++] = stress[i][1];
    buf[m++] = stress[i][2];
    buf[m++] = stress[i][3];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeActivestressAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    stress[j][0] += buf[m++];
    stress[j][1] += buf[m++];
    stress[j][2] += buf[m++];
    stress[j][3] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeActivestressAtom::memory_usage()
{
  double bytes = nmax*4 * sizeof(double);
  return bytes;
}
