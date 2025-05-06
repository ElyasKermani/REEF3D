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
Author: [Your Name]
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

// A specialized grid for efficient collision detection
class sixdof_collision_grid
{
public:
    sixdof_collision_grid(lexer *p, ghostcell *pgc);
    ~sixdof_collision_grid();
    
    // Update object positions in the grid
    void update_grid(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj);
    
    // Find potential collision pairs (broad phase)
    std::vector<std::pair<int, int>> find_potential_collisions(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj);
    
private:
    // Cell size for the grid (typically largest object diameter)
    double cell_size;
    
    // Domain boundaries
    double x_min, x_max, y_min, y_max, z_min, z_max;
    
    // Grid dimensions
    int nx, ny, nz;
    
    // Hash function to convert 3D index to 1D
    std::size_t hash_position(int i, int j, int k) const;
    
    // Calculate cell indices for an object
    std::tuple<int, int, int> calculate_cell_indices(const Eigen::Vector3d &position) const;
    
    // Map from cell hash to objects in that cell
    std::unordered_map<std::size_t, std::vector<int>> grid_cells;
};

#endif 