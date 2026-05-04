/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef CONTACT_HISTORY_H_
#define CONTACT_HISTORY_H_

#include<Eigen/Dense>

// Persistent per-pair tangential history used by contact-force kernels
struct ContactHistory
{
    Eigen::Vector3d s_t;   // tangential slip/displacement history at contact (m)
    bool in_contact;       // true while the pair is overlapping this step
    double t_last;         // simtime (s) when this history entry was last advanced
};

#endif
