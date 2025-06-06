/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sixdof_cfd::isource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_cfd::jsource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_cfd::ksource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_cfd::isource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_cfd::jsource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_cfd::ksource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_cfd::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void sixdof_cfd::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}
