#include<vector>
#include<array>
#include<cmath>
#include<iostream>
#include<fstream>
#include<matrix.h>

// mesh resolution
const int ni = 2; // number of cells in x-direction
const int nj = 2;
const int nk = 1;

Vector3 p[ni+1][nj+1][nk+1];

struct Line
{
    Vector3 p0;
    Vector3 p1;

    Line() = default;

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

    IJK() = default;

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

    return i + j * (ni+1) + k * (ni+1) * (nj+1) + 0;
}

struct Surface
{
    std::array<IJK, 4> ijk;
    int bc;

    Surface(IJK ijk0, IJK ijk1, IJK ijk2, IJK ijk3)
    {
        ijk[0] = ijk0;
        ijk[1] = ijk1;
        ijk[2] = ijk2;
        ijk[3] = ijk3;
    }
};

struct Face
{
    std::array<IJK, 4> ijk;

    Face(IJK ijk0, IJK ijk1, IJK ijk2, IJK ijk3)
    {
        ijk[0] = ijk0;
        ijk[1] = ijk1;
        ijk[2] = ijk2;
        ijk[3] = ijk3;
    }
};

struct Cell
{
    std::array<IJK, 8> ijk;

    Cell() = default;

    Cell(IJK ijk0, IJK ijk1, IJK ijk2, IJK ijk3, IJK ijk4, IJK ijk5, IJK ijk6, IJK ijk7)
    {
        ijk[0] = ijk0;
        ijk[1] = ijk1;
        ijk[2] = ijk2;
        ijk[3] = ijk3;
        ijk[4] = ijk4;
        ijk[5] = ijk5;
        ijk[6] = ijk6;
        ijk[7] = ijk7;
    }
};

Cell cells[ni][nj][nk];

Surface surfaces[ni][nj][2] IJ_surfaces;
Surface surfaces[2][nj][nk] JK_surfaces;
Surface surfaces[ni][2][nk] IK_surfaces;

Face faces[ni][nj][nk-1] IJ_faces;
Face faces[ni-1][nj][nk] JK_faces;
Face faces[ni][nj-1][nk] IK_faces;

void make_cells()
{
    for (int k=0; k<nk; ++k)
    {
        for (int j=0; j<nj; ++j)
        {
            for (int i=0; i<ni; ++i)
            {
                cells[i][j][k] = Cell(IJK(i,j,k), IJK(i+1, j, k), IJK(i+1, j+1, k), IJK(i, j+1, k), IJK(i,j,k+1), IJK(i+1, j, k+1), IJK(i+1, j+1, k+1), IJK(i, j+1, k+1));
            }
        }
    }
}

void make_surfaces()
{
    // IJ surfaces
    for (int j=0; j<nj; ++j)
    {
        for (int i=0; i<ni; ++i)
        {
            IJ_surfaces[i][j][0] = Surface(IJK(i,j,0), IJK(i+1, j, 0), IJK(i+1, j+1, 0), IJK(i, j+1, 0));
            IJ_surfaces[i][j][1] = Surface(IJK(i,j,nk), IJK(i+1, j, nk), IJK(i+1, j+1, nk), IJK(i, j+1, nk));
        }
    }

    // JK surfaces
    for (int j=0; j<nj; ++j)
    {
        for (int k=0; k<nk; ++k)
        {
            JK_surfaces[0][j][k] = Surface(IJK(0,j,k), IJK(0, j, k+1), IJK(0, j+1, k+1), IJK(0, j+1, k));
            JK_surfaces[1][j][k] = Surface(IJK(ni,j,k), IJK(ni, j, k+1), IJK(ni, j+1, k+1), IJK(ni, j+1, k));
        }
    }

    // IK surfaces
    for (int i=0; i<ni; ++i)
    {
        for (int k=0; k<nk; ++k)
        {
            IK_surfaces[i][0][k] = Surface(IJK(i,0,k), IJK(i+1, 0, k), IJK(i+1, 0, k+1), IJK(i, 0, k+1));
            IK_surfaces[i][1][k] = Surface(IJK(i,nj,k), IJK(i+1, nj, k), IJK(i+1, nj, k+1), IJK(i, nj, k+1));
        }
    }
}

void make_faces()
{
    // IJ faces
    for (int k=1; k<nk; ++k)
    {
        for (int j=0; j<nj; ++j)
        {
            for (int i=0; i<ni; ++i)
            {
                IJ_faces[i][j][k-1] = Face(IJK(i,j,k), IJK(i+1, j, k), IJK(i+1, j+1, k), IJK(i, j+1, k));
                faces.push_back(Face(IJK(i,j,k), IJK(i+1, j, k), IJK(i+1, j+1, k), IJK(i, j+1, k)));
                // parent: 
                // if k-1 > 0
                // 1: cell IJK(i,j,k-1)
                // 2: cell IJK(i,j,k)
            }
        }
    }

    // JK faces
    for (int i=1; i<ni; ++i)
    {
        for (int j=0; j<nj; ++j)
        {
            for (int k=0; k<nk; ++k)
            {
                faces.push_back(Face(IJK(i,j,k),  IJK(i, j, k+1),  IJK(i, j+1, k+1),  IJK(i, j+1, k)));

                // parent: 
                // if i-1 > 0
                // 1: cell IJK(i-1,j,k)
                // 2: cell IJK(i,j,k)
            }
        }
    }

    // IK faces
    for (int j=1; j<nj; ++j)
    {
        for (int i=0; i<ni; ++i)
        {
            for (int k=0; k<nk; ++k)
            {
                faces.push_back(Face(IJK(i,j,k),  IJK(i+1, j, k),  IJK(i+1, j, k+1),  IJK(i, j, k+1)));

                // parent: 
                // if j-1 > 0
                // 1: cell IJK(i,j-1,k)
                // 2: cell IJK(i,j,k)
            }
        }
    }
}

void generate_edge_points(Line* l)
{
    // discretize l0 and l2 and l4 and l6
    {
        auto pts0 = l[0].generate(ni);
        auto pts2 = l[2].generate(ni);
        auto pts4 = l[4].generate(ni);
        auto pts6 = l[6].generate(ni);
        for (int i=0; i<ni+1; ++i)
        {
            p[i][0][0] = pts0[i];
            p[i][nj][0] = pts2[i];
            p[i][0][nk] = pts4[i];
            p[i][nj][nk] = pts6[i];
        }
    }

    // discretize l1 and l3 and l5 and l7
    {
        auto pts1 = l[1].generate(nj);
        auto pts3 = l[3].generate(nj);
        auto pts5 = l[5].generate(nj);
        auto pts7 = l[7].generate(nj);
        for (int j=0; j<nj+1; ++j)
        {
            p[0][j][0] = pts3[j];
            p[ni][j][0] = pts1[j];
            p[0][j][nk] = pts5[j];
            p[ni][j][nk] = pts7[j];
        }
    }
}

void generate_surface_points()
{
    for (int i=1; i<ni; ++i)
    {
        for (int k=1; k<nk; ++k)
        {
            Line l(p[i][0][0], p[i][0][nk]);
            auto pts = l.generate(nk);
            p[i][0][k] = pts[k];
        }
    }
}

void generate_inner_points()
{
    for (int i=0; i<ni+1; ++i)
    {
        for (int k=0; k<nk+1; ++k)
        {
            Line l(p[i][0][k], p[i][nj][k]);
            auto pts = l.generate(nj);
            for (int j=1; j<nj; ++j)
            {
                p[i][j][k] = pts[j];
            }
        }
    }
}

void print_custom()
{
    int ncell = cells.size();
    int ncell_ij = 2 * ni * nj;
    int ncell_ik = 2 * ni * nk;
    int ncell_jk = 2 * nj * nk;
    int ncell_sur = ncell_ij + ncell_ik + ncell_jk;
    int ncell_all = ncell + ncell_ij + ncell_ik + ncell_jk;

    std::ofstream out;
    out.open("custom.mesh");

    out << "POINTS ";
    out << (ni+1) * (nj+1) * (nk+1) << std::endl;
    for (int k=0; k<nk+1; ++k)
    {
        for (int j=0; j<nj+1; ++j)
        {
            for (int i=0; i<ni+1; ++i)
            {
                out << p[i][j][k](0) << " " << p[i][j][k](1) << " " << p[i][j][k](2) << std::endl;
            }
        }
    }

    out << "CELLS ";
    out << cells.size() << std::endl;
    for (int k=0; k<nk; ++k)
    {
        for (int j=0; j<nj; ++j)
        {
            for (int i=0; i<ni; ++i)
            {
                out << cell.ijk[0] << " ";
                out << cell.ijk[1] << " ";
                out << cell.ijk[2] << " ";
                out << cell.ijk[3] << " ";
                out << cell.ijk[4] << " ";
                out << cell.ijk[5] << " ";
                out << cell.ijk[6] << " ";
                out << cell.ijk[7] << std::endl;
            }
        }
    }

    out << "IJ0 SURFACES ";
    out << IJ0_surfaces.size() << std::endl;

    out.close();
}

void print_vtk_str()
{
    int ncell = cells.size();
    int ncell_ij = 2 * ni * nj;
    int ncell_ik = 2 * ni * nk;
    int ncell_jk = 2 * nj * nk;
    int ncell_sur = ncell_ij + ncell_ik + ncell_jk;
    int ncell_all = ncell + ncell_ij + ncell_ik + ncell_jk;

    std::ofstream out;
    out.open("str.vtk");

    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "All in VTK format" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET STRUCTURED_GRID" << std::endl;
    out << "DIMENSIONS " << std::endl;
    out << ni+1;
    out << " ";
    out << nj+1;
    out << " ";
    out << nk+1 << std::endl;
    out << "POINTS ";
    out << (ni+1) * (nj+1) * (nk+1);
    out << " ";
    out << "float" << std::endl;
    for (int k=0; k<nk+1; ++k)
    {
        for (int j=0; j<nj+1; ++j)
        {
            for (int i=0; i<ni+1; ++i)
            {
                out << p[i][j][k](0) << " " << p[i][j][k](1) << " " << p[i][j][k](2) << std::endl;
            }
        }
    }

    out.close();
}

void print_vtk_ustr()
{
    int ncell = cells.size();
    int ncell_ij = 2 * ni * nj;
    int ncell_ik = 2 * ni * nk;
    int ncell_jk = 2 * nj * nk;
    int ncell_sur = ncell_ij + ncell_ik + ncell_jk;
    int ncell_all = ncell + ncell_ij + ncell_ik + ncell_jk;

    std::ofstream out;
    out.open("ustr.vtk");

    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "All in VTK format" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;
    out << "POINTS ";
    out << (ni+1) * (nj+1) * (nk+1);
    out << " ";
    out << "float" << std::endl;
    for (int k=0; k<nk+1; ++k)
    {
        for (int j=0; j<nj+1; ++j)
        {
            for (int i=0; i<ni+1; ++i)
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
}

int main()
{
    // Corner points
    Vector3 p0(0.0, 0.0, 0.0);
    Vector3 p1(1.0, 0.0, 0.0);
    Vector3 p2(1.0, 1.0, 0.0);
    Vector3 p3(0.0, 1.0, 0.0);

    Vector3 p4(0.0, 0.0, 1.0);
    Vector3 p5(1.0, 0.0, 1.0);
    Vector3 p6(1.0, 1.0, 1.0);
    Vector3 p7(0.0, 1.0, 1.0);

    // generate edges
    Line* l = new Line[8];

    l[0] = Line(p0, p1);
    l[1] = Line(p1, p2);
    l[2] = Line(p3, p2);
    l[3] = Line(p0, p3);

    l[4] = Line(p4, p5);
    l[5] = Line(p5, p6);
    l[6] = Line(p7, p6);
    l[7] = Line(p4, p7);

    generate_edge_points(l);
    generate_surface_points();
    generate_inner_points();
    make_surfaces();
    make_cells();
    make_faces();
    print_vtk_str();
    //print_vtk_ustr();




    return 0;
}
