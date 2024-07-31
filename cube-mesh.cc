#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;

void createCube(Triangulation<3>& ventricle)
{
    constexpr int dim = 3; // Define the dimension as a constant

    // Define vertices of the cube
    const std::vector<Point<dim>> vertices = {
        Point<dim>(0.0, 0.0, 0.0), // Vertex 0
        Point<dim>(1.0, 0.0, 0.0), // Vertex 1
        Point<dim>(0.0, 1.0, 0.0), // Vertex 2
        Point<dim>(1.0, 1.0, 0.0), // Vertex 3
        Point<dim>(0.0, 0.0, 1.0), // Vertex 4
        Point<dim>(1.0, 0.0, 1.0), // Vertex 5
        Point<dim>(0.0, 1.0, 1.0), // Vertex 6        
        Point<dim>(1.0, 1.0, 1.0), // Vertex 7
        Point<dim>(0.0, 0.0, 2.0), // Vertex 8
        Point<dim>(1.0, 0.0, 2.0), // Vertex 9
        Point<dim>(0.0, 1.0, 2.0), // Vertex 10        
        Point<dim>(1.0, 1.0, 2.0)  // Vertex 11
    };

    // Define the cell (hexahedron) connecting the vertices
    const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>> cell_vertices = 
    {
        {{0, 1, 2, 3, 4, 5, 6, 7}},
        {{4, 5, 6, 7, 8, 9, 10, 11}}
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
}

void convertToTetrahedra(const Triangulation<3>& hex_mesh, Triangulation<3>& tet_mesh)
{
    constexpr int dim = 3;

    // Get the vertices from the original hexahedral mesh
    std::vector<Point<dim>> vertices = hex_mesh.get_vertices();

    // Create a list of tetrahedrons
    std::vector<std::array<int, 4>> cell_vertices;

    // Iterate over each cell (hexahedron) in the original mesh
    for (const auto &cell : hex_mesh.active_cell_iterators())
    {
        // Get the indices of the current hexahedron's vertices
        std::vector<unsigned int> hex_vertex_indices(GeometryInfo<dim>::vertices_per_cell);
        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
            hex_vertex_indices[i] = cell->vertex_index(i);
        }

        // Define the tetrahedra splitting for the hexahedron
        const std::vector<std::array<int, 4>> tetrahedra = {
            {{3,7,6,4}}, // Tetrahedron 1
            {{2,3,6,4}}, // Tetrahedron 2
            {{2,0,3,4}}, // Tetrahedron 3
            {{3,5,7,4}}, // Tetrahedron 4
            {{3,1,5,4}}, // Tetrahedron 5
            {{0,1,3,4}}  // Tetrahedron 6
        };

        // Print the tetrahedra elements
        std::cout << "Tetrahedra elements for current hexahedron:" << std::endl;
        for (const auto &tet : tetrahedra)
        {
            std::cout << "Tetrahedron with vertex indices: ";
            for (int idx : tet)
            {
                std::cout << hex_vertex_indices[idx] << " ";
            }
            std::cout << std::endl;

            cell_vertices.push_back({{
                hex_vertex_indices[tet[0]],
                hex_vertex_indices[tet[1]],
                hex_vertex_indices[tet[2]],
                hex_vertex_indices[tet[3]]
            }});
        }
    }

    const unsigned int n_cells = cell_vertices.size();
    std::vector<CellData<dim>> cells(n_cells, CellData<dim>(4));
    for (unsigned int i = 0; i < n_cells; ++i)
    {
        for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
            cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
    }

    std::cout << "tet_cells.size() = " << cell_vertices.size() << std::endl;
    // Create the new triangulation
    tet_mesh.create_triangulation(vertices, cells, SubCellData());
}

int main()
{
    Triangulation<3> hex_mesh;
    createCube(hex_mesh);

    Triangulation<3> tet_mesh;
    convertToTetrahedra(hex_mesh, tet_mesh);

    std::ofstream out("tet_mesh.vtk");
    GridOut grid_out;
    grid_out.write_vtk(tet_mesh, out);

    std::cout << "Tetrahedral mesh has been written to tet_mesh.vtk" << std::endl;

    return 0;
}
