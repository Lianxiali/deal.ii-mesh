#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;

int main()
{
    constexpr int dim = 3; // Define the dimension as a constant

    Triangulation<dim> ventricle;

    // Define vertices of the cube
    const std::vector<Point<dim>> vertices = {
        Point<dim>(0.0, 0.0, 0.0), // Vertex 0
        Point<dim>(1.0, 0.0, 0.0), // Vertex 1
        Point<dim>(1.0, 1.0, 0.0), // Vertex 2
        Point<dim>(0.0, 1.0, 0.0), // Vertex 3
        Point<dim>(0.0, 0.0, 1.0), // Vertex 4
        Point<dim>(1.0, 0.0, 1.0), // Vertex 5
        Point<dim>(1.0, 1.0, 1.0), // Vertex 6
        Point<dim>(0.0, 1.0, 1.0),  // Vertex 7
        Point<dim>(0.0, 0.0, 2.0), // Vertex 8
        Point<dim>(1.0, 0.0, 2.0), // Vertex 9
        Point<dim>(1.0, 1.0, 2.0), // Vertex 10
        Point<dim>(0.0, 1.0, 2.0)  // Vertex 11
    };

    // Define the cell (hexahedron) connecting the vertices
    const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>>
    cell_vertices = 
        {
          {{0, 1, 3, 2, 4, 5, 7, 6}},
          {{4, 5, 7, 6, 8, 9, 11, 10}}
        };
    const unsigned int n_cells = cell_vertices.size();

    std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

    // Create the triangulation
    ventricle.create_triangulation(vertices, cells, SubCellData());

    // Output the grid to a file for visualization
    std::ofstream out("cube.vtk");
    GridOut grid_out;
    grid_out.write_vtk(ventricle, out);

    std::cout << "Cube mesh has been written to cube.vtk" << std::endl;

    return 0;
}
