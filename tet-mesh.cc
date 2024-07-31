#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;
/*
       6--------7
      /|       /|
     4--------5 |
     | |      | |
     | 2------|-3
     |/       |/
     0--------1
*/
void create_tet_mesh(Triangulation<3>& tet_mesh)
{
    constexpr int dim = 3;

    const std::vector<Point<dim>> vertices = {
        Point<dim>(0.0, 0.0, 0.0), // Vertex 0
        Point<dim>(1.0, 0.0, 0.0), // Vertex 1
        Point<dim>(0.0, 1.0, 0.0), // Vertex 2
        Point<dim>(1.0, 1.0, 0.0), // Vertex 3

        Point<dim>(0.0, 0.0, 1.0), // Vertex 4
        Point<dim>(1.0, 0.0, 1.0), // Vertex 5
        Point<dim>(0.0, 1.0, 1.0), // Vertex 6        
        Point<dim>(1.0, 1.0, 1.0)  // Vertex 7
    };

    // Define the cell (hexahedron) connecting the vertices
    const std::vector<std::array<int, 4>> cell_vertices = 
    {
            {{3,7,6,4}},
            {{2,3,6,4}},
            {{2,0,3,4}},
            {{3,5,7,4}},
            {{3,1,5,4}},
            {{0,1,3,4}}
    };

    const unsigned int n_cells = cell_vertices.size();

    std::vector<CellData<dim>> cells(n_cells, CellData<dim>(4));
    for (unsigned int i = 0; i < n_cells; ++i)
    {
        std::cout << "cellID = " << i << " : ";
        for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
        {
            cells[i].vertices[j] = cell_vertices[i][j];
            std::cout << cells[i].vertices[j];
        }
        std::cout << std::endl;
        cells[i].material_id = 0;
    }
std::cout << "cells.size() = " << cells.size() << std::endl;
    // Create the new triangulation
    tet_mesh.create_triangulation(vertices, cells, SubCellData());
}

int main()
{
    Triangulation<3> tet_mesh;
    create_tet_mesh(tet_mesh);

    std::ofstream out("tet_mesh.vtk");
    GridOut grid_out;
    grid_out.write_vtk(tet_mesh, out);

    std::cout << "Tetrahedral mesh has been written to tet_mesh.vtk" << std::endl;

    return 0;
}
