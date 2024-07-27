#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;

std::vector<Point<3>> createPointsOfCircle(int n, double radius, double z)
{
    constexpr int dim = 3; // Define the dimension as 3 (for a circle)

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

template<int dim>
void createVentricle
(
  const double& endo_short_radius, // outer
  const double& endo_long_radius, // inner
  const double& thickness,
  const double& base_elevation,
  const int& n_circumferential_cells,
  const int& n_radial_cells,
  const int& n_axial_cells,
  Triangulation<dim>& ventricle // Pass by reference to modify
)
{
    double epi_short_radius = endo_short_radius + thickness;
    double epi_long_radius  = endo_long_radius  + thickness;
    
    const int num_vertices_x = n_circumferential_cells;
    const int num_vertices_y = n_radial_cells + 1;
    const int num_vertices_z = n_axial_cells + 1;
    std::vector<Point<dim>> vertices;
    double u_max0 = -acos(base_elevation/endo_long_radius);
    double u_max1 = -acos(base_elevation/epi_long_radius);

    double u_min = -M_PI; // x0.999 to leave a tiny hole

    double du0 = (u_max0 - u_min) / n_axial_cells;
    double du1 = (u_max1 - u_min) / n_axial_cells;

    for(int k = 0; k < num_vertices_z; ++k)
    {      
      double u0 = u_min + k*du0;
      double u1 = u_min + k*du1;

      double z0 = endo_long_radius * cos(u0);
      double z1 = epi_long_radius  * cos(u1);

      double radius_inner = fabs(endo_short_radius*sin(u0)); // Radius of the circle
      double radius_outer = fabs(epi_short_radius*sin(u1));
      std::cout << "z = " << z0 << " z1 = " << z1 << std::endl;

      double dr = (radius_outer - radius_inner) / n_radial_cells;
      double dz = (z1 - z0)/n_radial_cells;

      for(int j = 0; j < num_vertices_y; ++j)
      {
        double radius = radius_inner + (j)*dr;
        double z = z0 + j * dz;

        std::cout << "radius = " << radius << std::endl;
        std::vector<Point<dim>> circle_points = createPointsOfCircle(num_vertices_x, radius, z);
        vertices.insert(vertices.end(), circle_points.begin(), circle_points.end());
      }
    }

    std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>> cell_vertices;
   
    for (int z = 0; z < n_axial_cells; ++z)
    {
        for (int y = 0; y < n_radial_cells; ++y)
        {
            for (int x = 0; x < n_circumferential_cells; ++x)
            {
                int v0 = x + y * num_vertices_x + z * num_vertices_x * num_vertices_y;
                int v1 = v0 + num_vertices_x;
                int v2 = v0 + 1;
                int v3 = v1 + 1;
                if(n_circumferential_cells-1 == x)// last closed cell shares some vertices of the first cell
                {
                  v2 = v0+1-n_circumferential_cells;
                  v3 = v1+1-n_circumferential_cells;
                }
                int v4 = v0 + num_vertices_x * num_vertices_y;        
                int v5 = v1 + num_vertices_x * num_vertices_y;
                int v6 = v2 + num_vertices_x * num_vertices_y;
                int v7 = v3 + num_vertices_x * num_vertices_y;

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

    ventricle.create_triangulation(vertices, cells, SubCellData());    
}

int main()
{
    double endo_short_radius = 7; // Radius of the circle
    double endo_long_radius = 17;
    const int n_circumferential_cells = 10;
    const int n_radial_cells = 3;
    const int n_axial_cells = 10; 
    const int thickness = 3;
    const double base_elevation = 5;

    constexpr int dim = 3; // Dimension of the space
    
    // Create the triangulation
    Triangulation<dim> ventricle;
    createVentricle<dim>(
        endo_short_radius,
        endo_long_radius,
        thickness,
        base_elevation,
        n_circumferential_cells,
        n_radial_cells,
        n_axial_cells,
        ventricle
    );

    // ventricle.refine_global(2);

    // Output the grid to a file for visualization
    std::ofstream out("ventricle.vtk");
    GridOut grid_out;
    grid_out.write_vtk(ventricle, out);

    std::cout << "Ventricle mesh has been written to ventricle.vtk" << std::endl;

    return 0;
}
