#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<matrix.h>

const int nx = 2; // number of cells in x-direction
const int ny = 2;
const int nz = 1;

Vector3 p[nx+1][ny+1][nz+1];

struct Line
{
    Vector3 p0;
    Vector3 p1;

    Line(Vector3& p0, Vector3& p1)
    {
        this->p0 = p0;
        this->p1 = p1;
    }

    double len()
    {
        Vector3 c = p1 - p0;
        return c.len();
    }

    double dl(int n)
    {
        return len() / double(n);
    }

    double dx(int n)
    {
        return std::abs(p1(0) - p0(0)) / double(n);
    }

    double dy(int n)
    {
        return std::abs(p1(1) - p0(1)) / double(n);
    }

    double dz(int n)
    {
        return std::abs(p1(2) - p0(2)) / double(n);
    }

    std::vector<Vector3> generate(int n)
    {
        std::vector<Vector3> pts;

        pts.push_back(p0);

        for (int i=0; i<n-1; ++i)
        {
            pts.push_back(pts.back() + Vector3(dx(n), dy(n), dz(n)));
        }

        pts.push_back(p1);

        return pts;
    }
};

struct IJK
{
    int i;
    int j;
    int k;

    IJK(int i, int j, int k)
    {
        this->i = i;
        this->j = j;
        this->k = k;
    }
};

int pos(IJK ijk)
{
    int i = ijk.i;
    int j = ijk.j;
    int k = ijk.k;

    return i + j * (nx+1) + k * (nx+1) * (ny+1) + 0;
}

struct Surface
{
    std::vector<IJK> ijk;

    Surface(std::vector<IJK> ijk)
    {
        this->ijk = ijk;
    }
};

struct Cell
{
    std::vector<IJK> ijk;

    Cell(std::vector<IJK> ijk)
    {
        this->ijk = ijk;
    }
};

std::vector<Cell> cells;
std::vector<Surface> surfaces;
std::vector<Surface> faces;

void make_cells()
{
    for (int k=0; k<nz; ++k)
    {
        for (int j=0; j<ny; ++j)
        {
            for (int i=0; i<nx; ++i)
            {
                cells.push_back(std::vector<IJK>{IJK(i,j,k), IJK(i+1, j, k), IJK(i+1, j+1, k), IJK(i, j+1, k), IJK(i,j,k+1), IJK(i+1, j, k+1), IJK(i+1, j+1, k+1), IJK(i, j+1, k+1)});
            }
        }
    }
}

void make_surfaces()
{
    for (int j=0; j<ny; ++j)
    {
        for (int i=0; i<nx; ++i)
        {
            surfaces.push_back(std::vector<IJK>{IJK(i,j,0),  IJK(i+1, j, 0),  IJK(i+1, j+1, 0),  IJK(i, j+1, 0)});
            surfaces.push_back(std::vector<IJK>{IJK(i,j,nz), IJK(i+1, j, nz), IJK(i+1, j+1, nz), IJK(i, j+1, nz)});
        }
    }

    for (int j=0; j<ny; ++j)
    {
        for (int k=0; k<nz; ++k)
        {
            surfaces.push_back(std::vector<IJK>{IJK(0,j,k),  IJK(0, j, k+1),  IJK(0, j+1, k+1),  IJK(0, j+1, k)});
            surfaces.push_back(std::vector<IJK>{IJK(nx,j,k),  IJK(nx, j, k+1),  IJK(nx, j+1, k+1),  IJK(nx, j+1, k)});
        }
    }

    for (int i=0; i<nx; ++i)
    {
        for (int k=0; k<nz; ++k)
        {
            surfaces.push_back(std::vector<IJK>{IJK(i,0,k),  IJK(i+1, 0, k),  IJK(i+1, 0, k+1),  IJK(i, 0, k+1)});
            surfaces.push_back(std::vector<IJK>{IJK(i,ny,k),  IJK(i+1, ny, k),  IJK(i+1, ny, k+1),  IJK(i, ny, k+1)});
        }
    }
}

void make_faces()
{
    for (int k=1; k<nz; ++k)
    {
        for (int j=0; j<ny; ++j)
        {
            for (int i=0; i<nx; ++i)
            {
                faces.push_back(std::vector<IJK>{IJK(i,j,k), IJK(i+1, j, k), IJK(i+1, j+1, k), IJK(i, j+1, k)});
		// parent: 
		// if k-1 > 0
		// 1: cell IJK(i,j,k-1)
		// 2: cell IJK(i,j,k)
            }
        }
    }

    for (int i=1; i<nx; ++i)
    {
	    for (int j=0; j<ny; ++j)
	    {
		    for (int k=0; k<nz; ++k)
		    {
			    faces.push_back(std::vector<IJK>{IJK(i,j,k),  IJK(i, j, k+1),  IJK(i, j+1, k+1),  IJK(i, j+1, k)});

		// parent: 
		// if i-1 > 0
		// 1: cell IJK(i-1,j,k)
		// 2: cell IJK(i,j,k)
		    }
	    }
    }

    for (int j=1; j<ny; ++j)
    {
	    for (int i=0; i<nx; ++i)
	    {
		    for (int k=0; k<nz; ++k)
		    {
			    faces.push_back(std::vector<IJK>{IJK(i,j,k),  IJK(i+1, j, k),  IJK(i+1, j, k+1),  IJK(i, j, k+1)});

		// parent: 
		// if j-1 > 0
		// 1: cell IJK(i,j-1,k)
		// 2: cell IJK(i,j,k)
		    }
	    }
    }

}

int main()
{
    Vector3 p0(0.0, 0.0, 0.0);
    Vector3 p1(1.0, 0.0, 0.0);
    Vector3 p2(1.0, 1.0, 0.0);
    Vector3 p3(0.0, 1.0, 0.0);

    Vector3 p4(0.0, 0.0, 1.0);
    Vector3 p5(1.0, 0.0, 1.0);
    Vector3 p6(1.0, 1.0, 1.0);
    Vector3 p7(0.0, 1.0, 1.0);

    Line l0(p0, p1);
    Line l1(p1, p2);
    Line l2(p3, p2);
    Line l3(p0, p3);

    Line l4(p4, p5);
    Line l5(p5, p6);
    Line l6(p7, p6);
    Line l7(p4, p7);

    Line l8(p0, p4);
    Line l9(p1, p5);
    Line l10(p2, p6);
    Line l11(p3, p7);


    // discretize l0 and l2 and l4 and l6
    {
        auto pts0 = l0.generate(nx);
        auto pts2 = l2.generate(nx);
        auto pts4 = l4.generate(nx);
        auto pts6 = l6.generate(nx);
        for (int i=0; i<nx+1; ++i)
        {
            p[i][0][0] = pts0[i];
            p[i][ny][0] = pts2[i];
            p[i][0][nz] = pts4[i];
            p[i][ny][nz] = pts6[i];
        }
    }

    // discretize l1 and l3 and l5 and l7
    {
        auto pts1 = l1.generate(ny);
        auto pts3 = l3.generate(ny);
        auto pts5 = l5.generate(ny);
        auto pts7 = l7.generate(ny);
        for (int j=0; j<ny+1; ++j)
        {
            p[0][j][0] = pts3[j];
            p[nx][j][0] = pts1[j];
            p[0][j][nz] = pts5[j];
            p[nx][j][nz] = pts7[j];
        }
    }

    // discretize l8 and l9 and l10 and l11
    {
        auto pts8 = l8.generate(nz);
        auto pts9 = l9.generate(nz);
        auto pts10 = l10.generate(nz);
        auto pts11 = l11.generate(nz);
        for (int k=0; k<nz+1; ++k)
        {
            p[0][0][k] = pts8[k];
            p[nx][0][k] = pts9[k];
            p[nx][ny][k] = pts10[k];
            p[0][ny][k] = pts11[k];
        }
    }

    {
        for (int i=1; i<nx; ++i)
        {
            for (int k=1; k<nz; ++k)
            {
                Line l(p[i][0][0], p[i][0][nz]);
                auto pts = l.generate(nz);
                p[i][0][k] = pts[k];
            }
        }
    }

    {
        for (int i=1; i<nx; ++i)
        {
            for (int k=1; k<nz; ++k)
            {
                Line l(p[i][ny][0], p[i][ny][nz]);
                auto pts = l.generate(nz);
                p[i][ny][k] = pts[k];
            }
        }
    }

    for (int i=0; i<nx+1; ++i)
    {
        for (int k=0; k<nz+1; ++k)
        {
            Line l(p[i][0][k], p[i][ny][k]);
            auto pts = l.generate(ny);
            for (int j=0; j<ny+1; ++j)
            {
                p[i][j][k] = pts[j];
            }
        }
    }

    make_surfaces();
    make_cells();
    make_faces();

    int ncell = cells.size();
    int ncell_xy = 2 * nx * ny;
    int ncell_xz = 2 * nx * nz;
    int ncell_yz = 2 * ny * nz;
    int ncell_sur = ncell_xy + ncell_xz + ncell_yz;
    int ncell_all = ncell + ncell_xy + ncell_xz + ncell_yz;

    std::cout << "ncell: " << ncell << std::endl;
    std::cout << "ncell_xy: " << ncell_xy << std::endl;
    std::cout << "ncell_xz: " << ncell_xz << std::endl;
    std::cout << "ncell_yz: " << ncell_yz << std::endl;
    
    std::ofstream out;
    out.open("str.vtk");

    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "All in VTK format" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET STRUCTURED_GRID" << std::endl;
    out << "DIMENSIONS " << std::endl;
    out << nx+1;
    out << " ";
    out << ny+1;
    out << " ";
    out << nz+1 << std::endl;
    out << "POINTS ";
    out << (nx+1) * (ny+1) * (nz+1);
    out << " ";
    out << "float" << std::endl;
    for (int k=0; k<nz+1; ++k)
    {
        for (int j=0; j<ny+1; ++j)
        {
            for (int i=0; i<nx+1; ++i)
            {
                out << p[i][j][k](0) << " " << p[i][j][k](1) << " " << p[i][j][k](2) << std::endl;
            }
        }
    }

    out.close();


    out.open("ustr.vtk");

    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "All in VTK format" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;
    out << "POINTS ";
    out << (nx+1) * (ny+1) * (nz+1);
    out << " ";
    out << "float" << std::endl;
    for (int k=0; k<nz+1; ++k)
    {
        for (int j=0; j<ny+1; ++j)
        {
            for (int i=0; i<nx+1; ++i)
            {
                out << p[i][j][k](0) << " " << p[i][j][k](1) << " " << p[i][j][k](2) << std::endl;
            }
        }
    }
    out << "CELLS ";
    out << ncell_all;
    out << " ";
    out << cells.size() * 9 + surfaces.size() * 5 << std::endl;
    for (Surface& surface: surfaces)
    {
        out << 4 << " ";
	out << pos(surface.ijk[0]) << " ";
	out << pos(surface.ijk[1]) << " ";
	out << pos(surface.ijk[2]) << " ";
	out << pos(surface.ijk[3]) << " ";
        out << std::endl;
    }

    for (Cell& cell: cells)
    {
        out << 8 << " ";
	out << pos(cell.ijk[0]) << " ";
	out << pos(cell.ijk[1]) << " ";
	out << pos(cell.ijk[2]) << " ";
	out << pos(cell.ijk[3]) << " ";
	out << pos(cell.ijk[4]) << " ";
	out << pos(cell.ijk[5]) << " ";
	out << pos(cell.ijk[6]) << " ";
	out << pos(cell.ijk[7]) << " ";
	out << std::endl;
    }

    out << "CELL_TYPES ";
    out << ncell_all << std::endl;
    for (Surface& surface: surfaces)
    {
        out << 9 << std::endl;
    }
    for (Cell& cell: cells)
    {
        out << 12 << std::endl;
    }

    out.close();

    return 0;
}
