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

#include"sixdof_collision_grid.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<iostream>
#include<algorithm>
#include<cmath>

sixdof_collision_grid::sixdof_collision_grid(lexer *p, ghostcell *pgc)
{
    // Set domain boundaries based on computational domain
    x_min = p->originx;
    x_max = p->originx + p->xmax*p->DXM;
    y_min = p->originy;
    y_max = p->originy + p->ymax*p->DYM;
    z_min = p->originz;
    z_max = p->originz + p->zmax*p->DZM;
    
    // Default cell size - this can be adjusted based on object sizes
    cell_size = 2.0 * p->DXM; // Start with twice the grid cell size
    
    // Calculate grid dimensions
    nx = static_cast<int>(std::ceil((x_max - x_min) / cell_size));
    ny = static_cast<int>(std::ceil((y_max - y_min) / cell_size));
    nz = static_cast<int>(std::ceil((z_max - z_min) / cell_size));
    
    if(p->mpirank==0)
    {
        std::cout<<"6DOF Collision Grid initialized..."<<std::endl;
        std::cout<<"  Grid dimensions: "<<nx<<" x "<<ny<<" x "<<nz<<std::endl;
        std::cout<<"  Cell size: "<<cell_size<<std::endl;
    }
}

sixdof_collision_grid::~sixdof_collision_grid()
{
}

void sixdof_collision_grid::update_grid(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj)
{
    // Clear the grid
    grid_cells.clear();
    
    // Insert each object into the grid
    for(int i=0; i<p->X20; ++i)
    {
        // Get object center
        Eigen::Vector3d center = fb_obj[i]->c_;
        
        // Calculate cell indices
        auto [ci, cj, ck] = calculate_cell_indices(center);
        
        // Calculate hash
        std::size_t hash = hash_position(ci, cj, ck);
        
        // Insert object index into the grid cell
        grid_cells[hash].push_back(i);
        
        // Also insert into neighboring cells based on object radius
        // This ensures objects that span multiple cells are properly detected
        double radius = fb_obj[i]->radius;
        
        // Calculate cell span based on radius
        int cell_span = static_cast<int>(std::ceil(radius / cell_size));
        
        // Insert into neighboring cells within the span
        for(int di = -cell_span; di <= cell_span; ++di)
        {
            for(int dj = -cell_span; dj <= cell_span; ++dj)
            {
                for(int dk = -cell_span; dk <= cell_span; ++dk)
                {
                    // Skip the center cell (already added)
                    if(di == 0 && dj == 0 && dk == 0)
                        continue;
                    
                    // Calculate neighboring cell indices
                    int ni = ci + di;
                    int nj = cj + dj;
                    int nk = ck + dk;
                    
                    // Skip if out of bounds
                    if(ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
                        continue;
                    
                    // Calculate hash for neighbor
                    std::size_t neighbor_hash = hash_position(ni, nj, nk);
                    
                    // Insert object index into the neighbor cell
                    grid_cells[neighbor_hash].push_back(i);
                }
            }
        }
    }
}

std::vector<std::pair<int, int>> sixdof_collision_grid::find_potential_collisions(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj)
{
    std::vector<std::pair<int, int>> collision_pairs;
    
    // Iterate through all populated grid cells
    for(const auto &cell_entry : grid_cells)
    {
        const auto &objects_in_cell = cell_entry.second;
        
        // If cell has more than one object, check for potential collisions
        if(objects_in_cell.size() > 1)
        {
            // Check all pairs of objects in this cell
            for(size_t i = 0; i < objects_in_cell.size() - 1; ++i)
            {
                for(size_t j = i + 1; j < objects_in_cell.size(); ++j)
                {
                    int obj1_index = objects_in_cell[i];
                    int obj2_index = objects_in_cell[j];
                    
                    // Avoid duplicate pairs (objects might be in multiple cells)
                    if(std::find(collision_pairs.begin(), collision_pairs.end(), 
                                std::make_pair(obj1_index, obj2_index)) == collision_pairs.end() &&
                       std::find(collision_pairs.begin(), collision_pairs.end(), 
                                std::make_pair(obj2_index, obj1_index)) == collision_pairs.end())
                    {
                        // Add this pair as a potential collision
                        collision_pairs.push_back(std::make_pair(obj1_index, obj2_index));
                    }
                }
            }
        }
    }
    
    return collision_pairs;
}

std::size_t sixdof_collision_grid::hash_position(int i, int j, int k) const
{
    // Simple hash function: i + j*nx + k*nx*ny
    return static_cast<std::size_t>(i + j * nx + k * nx * ny);
}

std::tuple<int, int, int> sixdof_collision_grid::calculate_cell_indices(const Eigen::Vector3d &position) const
{
    // Calculate grid cell indices from position
    int i = static_cast<int>((position(0) - x_min) / cell_size);
    int j = static_cast<int>((position(1) - y_min) / cell_size);
    int k = static_cast<int>((position(2) - z_min) / cell_size);
    
    // Clamp to valid range
    i = std::max(0, std::min(i, nx - 1));
    j = std::max(0, std::min(j, ny - 1));
    k = std::max(0, std::min(k, nz - 1));
    
    return {i, j, k};
} 