#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;

using namespace dealii;

std::vector<Point<3>> createPointsOfCircle(int n, double radius, double z)
{
    constexpr int dim = 3; // Define the dimension as 2 (for a circle)

    // Define angles for the n points
    std::vector<double> angles;
    double dsita = 2 * M_PI / n;
    for (int i = 0; i < n; ++i)
    {
        angles.push_back(i * dsita);
    }

    // Define points for the circle
    std::vector<Point<dim>> circle_points;
    for (const auto &angle : angles)
    {
        circle_points.push_back(Point<dim>(radius * std::cos(angle),
                                           radius * std::sin(angle),
                                           z));
    }

    return circle_points;    
}

int main()
{

    double rs0 = 7; // Radius of the circle
    double rs1 = 10;
    double rl0 = 17;
    double rl1 = 20;   
    const int num_cells_x = 40;
    const int num_cells_y = 3;
    const int num_cells_z = 50; 
    const int num_vertices_x = num_cells_x;
    const int num_vertices_y = num_cells_y+1;
    const int num_vertices_z = num_cells_z+1;

    constexpr int dim = 3; // Dimension of the space
    
    std::vector<Point<3>> vertices;
    double u_max0 = -acos(5.0/17.0);
    double u_max1 = -acos(5.0/20.0);

    double u_min = -M_PI;

    double du0 = (u_max0 - u_min) / num_cells_z;
    double du1 = (u_max1 - u_min) / num_cells_z;

    for(int k = 0; k < num_vertices_z; ++k)
    {
      
      double u0 = u_min + k*du0;
      double u1 = u_min + k*du1;

      double z0 = rl0 * cos(u0);
      double z1 = rl1 * cos(u1);

      double radius_inner = fabs(rs0*sin(u0)); // Radius of the circle
      double radius_outer = fabs(rs1*sin(u1));
      std::cout << "z = " << z0 << " z1 = " << z1 << std::endl;

      double dr = (radius_outer - radius_inner) / num_cells_y;
      double dz = (z1 - z0)/num_cells_y;

      for(int j = 0; j < num_vertices_y; ++j)
      {
        double radius = radius_inner + (j)*dr;
        double z = z0 + j * dz;

        std::cout << "radius = " << radius << std::endl;
        std::vector<Point<3>> circle_points = createPointsOfCircle(num_vertices_x, radius, z);
        vertices.insert(vertices.end(), circle_points.begin(), circle_points.end());
      }
    }

    std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>> cell_vertices;
   
    for (int z = 0; z < num_cells_z; ++z)
    {
        for (int y = 0; y < num_cells_y; ++y)
        {
            int num_vertex_per_layer = (num_vertices_x + num_vertices_x) * num_vertices_y;
            for (int x = 0; x < num_cells_x; ++x)
            {
                int v0 = x + y * num_vertices_x + z * num_vertices_x * num_vertices_y;
                int v1 = v0 + num_vertices_x;
                int v2 = v0 + 1;
                int v3 = v1 + 1;
                if(num_cells_x-1 == x)// last closed cell shares some vertices of the first cell
                {
                  v2 = v0+1-num_cells_x;
                  v3 = v1+1-num_cells_x;
                }
                int v4 = v0 + num_vertices_x * num_vertices_y;        
                int v5 = v1 + num_vertices_x * num_vertices_y;
                int v6 = v2 + num_vertices_x * num_vertices_y;
                int v7 = v3 + num_vertices_x * num_vertices_y;

                // std::cout << v0 << " "
                //           << v1 << " "
                //           << v2 << " "
                //           << v3 << " "
                //           << v4 << " "
                //           << v5 << " "
                //           << v6 << " "
                //           << v7 << std::endl;

                cell_vertices.push_back({{v0, v1, v2, v3, v4, v5, v6, v7}});
            }
        }
    }

    const unsigned int n_cells = cell_vertices.size();

    std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

    // Create the triangulation
    Triangulation<dim> ventricle;
    ventricle.create_triangulation(vertices, cells, SubCellData());

    // Output the grid to a file for visualization
    std::ofstream out("ventricle.vtk");
    GridOut grid_out;
    grid_out.write_vtk(ventricle, out);

    std::cout << "Ventricle mesh has been written to ventricle.vtk" << std::endl;

    return 0;
}
