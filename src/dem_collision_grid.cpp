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

#include"dem_collision_grid.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<iostream>
#include<algorithm>
#include<cmath>

dem_collision_grid::dem_collision_grid(lexer *p, ghostcell *pgc)
{
    x_min = p->originx;
    x_max = p->originx + p->xmax*p->DXM;
    y_min = p->originy;
    y_max = p->originy + p->ymax*p->DYM;
    z_min = p->originz;
    z_max = p->originz + p->zmax*p->DZM;

    // Initial cell size; the first update_grid() call adapts it to the largest object
    cell_size = 2.0 * p->DXM;

    nx = static_cast<int>(std::ceil((x_max - x_min) / cell_size));
    ny = static_cast<int>(std::ceil((y_max - y_min) / cell_size));
    nz = static_cast<int>(std::ceil((z_max - z_min) / cell_size));

    if(p->mpirank==0)
    {
        std::cout<<"6DOF Collision Grid initialized..."<<std::endl;
        std::cout<<"  Grid dimensions: "<<nx<<" x "<<ny<<" x "<<nz<<std::endl;
        std::cout<<"  Initial cell size: "<<cell_size<<std::endl;
        std::cout<<"  Note: Cell size will be adapted based on object radii"<<std::endl;
    }
}

dem_collision_grid::~dem_collision_grid()
{
}

double dem_collision_grid::compute_optimal_cell_size(lexer *p, std::vector<sixdof_obj*> &fb_obj)
{
    double max_radius = 0.0;
    for(int i = 0; i < p->X20; ++i)
    {
        if(fb_obj[i]->radius > max_radius)
            max_radius = fb_obj[i]->radius;
    }

    // Aim for cells ~2.5 times the largest body so most objects fit in one cell,
    // but never go below twice the underlying fluid grid spacing.
    double optimal_size = 2.5 * max_radius;
    double min_size     = 2.0 * std::max({p->DXM, p->DYM, p->DZM});

    return std::max(optimal_size, min_size);
}

void dem_collision_grid::update_grid(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj)
{
    // Re-pick the cell size only when the optimum has drifted by more than 10%
    double new_cell_size = compute_optimal_cell_size(p, fb_obj);
    if(std::abs(new_cell_size - cell_size) > 0.1 * cell_size)
    {
        cell_size = new_cell_size;
        nx = static_cast<int>(std::ceil((x_max - x_min) / cell_size));
        ny = static_cast<int>(std::ceil((y_max - y_min) / cell_size));
        nz = static_cast<int>(std::ceil((z_max - z_min) / cell_size));

        if(p->mpirank == 0 && p->count % p->P12 == 0)
        {
            std::cout << "6DOF Collision Grid: Adapted cell size to " << cell_size << std::endl;
            std::cout << "  New grid dimensions: " << nx << " x " << ny << " x " << nz << std::endl;
        }
    }

    grid_cells.clear();

    // Insert each object into its centre cell and all neighbours covered by its radius
    for(int i=0; i<p->X20; ++i)
    {
        Eigen::Vector3d center = fb_obj[i]->c_;

        auto [ci, cj, ck] = calculate_cell_indices(center);
        grid_cells[hash_position(ci, cj, ck)].push_back(i);

        int cell_span = static_cast<int>(std::ceil(fb_obj[i]->radius / cell_size));

        for(int di = -cell_span; di <= cell_span; ++di)
        {
            for(int dj = -cell_span; dj <= cell_span; ++dj)
            {
                for(int dk = -cell_span; dk <= cell_span; ++dk)
                {
                    if(di == 0 && dj == 0 && dk == 0)
                        continue;

                    int ni = ci + di;
                    int nj = cj + dj;
                    int nk = ck + dk;

                    if(ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
                        continue;

                    grid_cells[hash_position(ni, nj, nk)].push_back(i);
                }
            }
        }
    }
}

std::vector<std::pair<int, int>> dem_collision_grid::find_potential_collisions(lexer *p, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj)
{
    std::vector<std::pair<int, int>> collision_pairs;

    for(const auto &cell_entry : grid_cells)
    {
        const auto &objects_in_cell = cell_entry.second;

        if(objects_in_cell.size() > 1)
        {
            for(size_t i = 0; i < objects_in_cell.size() - 1; ++i)
            {
                for(size_t j = i + 1; j < objects_in_cell.size(); ++j)
                {
                    int obj1_index = objects_in_cell[i];
                    int obj2_index = objects_in_cell[j];

                    // Objects can appear in several cells when their radius crosses
                    // a boundary; deduplicate the pair list before returning it.
                    if(std::find(collision_pairs.begin(), collision_pairs.end(),
                                std::make_pair(obj1_index, obj2_index)) == collision_pairs.end() &&
                       std::find(collision_pairs.begin(), collision_pairs.end(),
                                std::make_pair(obj2_index, obj1_index)) == collision_pairs.end())
                    {
                        collision_pairs.push_back(std::make_pair(obj1_index, obj2_index));
                    }
                }
            }
        }
    }

    return collision_pairs;
}

std::size_t dem_collision_grid::hash_position(int i, int j, int k) const
{
    return static_cast<std::size_t>(i + j * nx + k * nx * ny);
}

std::tuple<int, int, int> dem_collision_grid::calculate_cell_indices(const Eigen::Vector3d &position) const
{
    int i = static_cast<int>((position(0) - x_min) / cell_size);
    int j = static_cast<int>((position(1) - y_min) / cell_size);
    int k = static_cast<int>((position(2) - z_min) / cell_size);

    // Clamp to a valid grid index for objects sitting on the boundary
    i = std::max(0, std::min(i, nx - 1));
    j = std::max(0, std::min(j, ny - 1));
    k = std::max(0, std::min(k, nz - 1));

    return {i, j, k};
}