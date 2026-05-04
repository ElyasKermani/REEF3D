/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef SIXDOF_COLLISION_GRID_H_
#define SIXDOF_COLLISION_GRID_H_

#include<vector>
#include<map>
#include<unordered_map>
#include<tuple>
#include<Eigen/Dense>

class lexer;
class ghostcell;
class sixdof_obj;

// Spatial-hash broad-phase used by the 6DOF collision pipeline
class sixdof_collision_grid
{
public:
    sixdof_collision_grid(lexer *p, ghostcell *pgc);
    ~sixdof_collision_grid();

    // Bin every object into the grid for the current step
    void update_grid(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj);

    // Return the unique pairs of objects that share a cell
    std::vector<std::pair<int, int>> find_potential_collisions(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj);

    // Pick a cell size large enough to contain the biggest object
    double compute_optimal_cell_size(lexer *p, std::vector<sixdof_obj*> &fb_obj);

private:
    double cell_size;
    double x_min, x_max, y_min, y_max, z_min, z_max;
    int nx, ny, nz;

    // Linearize a 3D cell index into a hash key
    std::size_t hash_position(int i, int j, int k) const;

    // Find the cell that contains a given world-space position
    std::tuple<int, int, int> calculate_cell_indices(const Eigen::Vector3d &position) const;

    // Cell hash -> list of object indices currently in that cell
    std::unordered_map<std::size_t, std::vector<int>> grid_cells;
};

#endif 