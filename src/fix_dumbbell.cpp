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
#include "stdio.h"
#include "string.h"
#include "fix_dumbbell.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "math.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CCW,CW,MIXED,CONVECT,COM};

/* ---------------------------------------------------------------------- */

// example command:
// fix dumbbell FORCE_MAGNITUDE FORCE_TYPE
FixDumbbell::FixDumbbell(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"dumbbell") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix dumbbell command: not enough args");
  f_active = force->numeric(FLERR,arg[3]);
  activestyle = CCW;
  if (narg == 5){
    if (strcmp(arg[4],"ccw") == 0)
      activestyle = CCW;
    else if (strcmp(arg[4],"cw") == 0)
      activestyle = CW;
    else if (strcmp(arg[4],"mixed") == 0)
      activestyle = MIXED;
    else if (strcmp(arg[4],"convect") == 0)
      activestyle = CONVECT;
    else if (strcmp(arg[4],"com") == 0)
      activestyle = COM;
    else
      error->all(FLERR, "Only {ccw, cw, mixed, convect, com} are accepted styles.");
  }
  if (force->newton_bond)
    error->all(FLERR, "To use fix dumbbell, you must turn off newton bonds "
               "in the input file, e.g. with 'newton on off'.");
}

/* ---------------------------------------------------------------------- */

int FixDumbbell::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDumbbell::post_force(int /*vflag*/)
{
  double delx, dely, rsq, r;
  double f1x, f1y, f2x, f2y;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;

  if (activestyle==COM){  // Rotate CCW a dumbbell that includes COM dummy atom
    for (int n = 0; n < nanglelist; n++) {
      int i1 = anglelist[n][0];
      int i2 = anglelist[n][1];
      // Get a unit vector pointing from atom 1 to atom 2 (assuming 2d in xy-plane)
      delx = x[i2][0] - x[i1][0];
      dely = x[i2][1] - x[i1][1];
      rsq = delx*delx + dely*dely;
      r = sqrt(rsq);
      delx /= r;
      dely /= r;
      // Apply forces for a net CCW torque.
      f1x = f_active * (dely);  // unit vector rotated CW
      f1y = f_active * (-delx);
      f2x = f_active * (-dely); // unit vector rotated CCW
      f2y = f_active * (delx);
      // Finally, add the computed forces to atoms owned by this processor:
      if (i1 < nlocal){   // Add force only to real atoms (not ghosts)
        f[i1][0] += f1x;  // unit vector rotated CW
        f[i1][1] += f1y;
      }
      if (i2 < nlocal){   // Add force only to real atoms (not ghosts)
        f[i2][0] += f2x;  // unit vector rotated CCW
        f[i2][1] += f2y;
      }
    }
  }
  else {

    // Add the active force to the already-computed per-atom forces
    for (int n = 0; n < nbondlist; n++) {
      int i1 = bondlist[n][0];
      int i2 = bondlist[n][1];
      // Get a unit vector pointing from atom 1 to atom 2 (assuming 2d in xy-plane)
      delx = x[i2][0] - x[i1][0];
      dely = x[i2][1] - x[i1][1];
      rsq = delx*delx + dely*dely;
      r = sqrt(rsq);
      delx /= r;
      dely /= r;

      if (activestyle==CCW){
        // Apply forces for a net CCW torque.
        f1x = f_active * (dely);  // unit vector rotated CW
        f1y = f_active * (-delx);
        f2x = f_active * (-dely); // unit vector rotated CCW
        f2y = f_active * (delx);
      }

      else if (activestyle==CW){
        // Apply forces for a net CW torque.
        f1x = f_active * (-dely); // unit vector rotated CCW
        f1y = f_active * (delx);
        f2x = f_active * (dely);  // unit vector rotated CW
        f2y = f_active * (-delx);
      }

      else if (activestyle==MIXED){
        if (n % 2 == 0){
        // Apply forces for a net CCW torque.
        f1x = f_active * (dely);  // unit vector rotated CW
        f1y = f_active * (-delx);
        f2x = f_active * (-dely); // unit vector rotated CCW
        f2y = f_active * (delx);
        }
        else{
        // Apply forces for a net CW torque.
        f1x = f_active * (-dely); // unit vector rotated CCW
        f1y = f_active * (delx);
        f2x = f_active * (dely);  // unit vector rotated CW
        f2y = f_active * (-delx);
        }
      }

      else if (activestyle==CONVECT){
        // Apply convective force along bond axis
        f1x = f_active * (delx);
        f1y = f_active * (dely);
        f2x = f_active * (delx);
        f2y = f_active * (dely);
      }

      // Finally, add the computed forces to atoms owned by this processor:
      if (i1 < nlocal){   // Add force only to real atoms (not ghosts)
        f[i1][0] += f1x;  // unit vector rotated CW
        f[i1][1] += f1y;
      }
      if (i2 < nlocal){   // Add force only to real atoms (not ghosts)
        f[i2][0] += f2x;  // unit vector rotated CCW
        f[i2][1] += f2y;
      }
    }
  }
}
