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
#include "compute_dumbbellangle.h"
#include "neighbor.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

// example command:
// compute angle all dumbbellangle/atom
ComputeDumbbellAngle::ComputeDumbbellAngle(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  angles(NULL)

{
  if (narg < 3) error->all(FLERR,"Illegal compute dumbbellangle/atom command");
  comflag = 0;
  if (narg == 4)
    {
      if (strcmp(arg[3],"com") == 0) comflag = 1;
    }
  peratom_flag = 1;
  size_peratom_cols = 0;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeDumbbellAngle::~ComputeDumbbellAngle()
{
  memory->destroy(angles);
}

/* ---------------------------------------------------------------------- */

void ComputeDumbbellAngle::compute_peratom()
{

  invoked_peratom = update->ntimestep;

  // grow angles array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(angles);
    nmax = atom->nmax;
    memory->create(angles,nmax,"dumbbellangle/atom:angles");
    vector_atom = angles;
  }


  double **x = atom->x;
  int i1, i2;
  double dx, dy, angle;
  int nlocal = atom->nlocal;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;

  if (comflag)
    {
      for (int n = 0; n < nanglelist; n++)
        {
          i1 = anglelist[n][0];
          i2 = anglelist[n][1];
          dx = x[i2][0] - x[i1][0];
          dy = x[i2][1] - x[i1][1];
          angle = atan(fabs(dy/dx));
          angles[i1] = angle;
          angles[i2] = angle;
        }
    }

  else
    {
      for (int n = 0; n < nbondlist; n++)
        {
          i1 = bondlist[n][0];
          i2 = bondlist[n][1];
          dx = x[i2][0] - x[i1][0];
          dy = x[i2][1] - x[i1][1];
          angle = atan(fabs(dy/dx));
          angles[i1] = angle;
          angles[i2] = angle;
        }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDumbbellAngle::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
