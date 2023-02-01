#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<matrix.h>

int nx = 2; // number of cells in x-direction
int ny = 2;
int nz = 1;

int pos(int i, int j, int k)
{
    return i + j * (nx+1) + k * (nx+1) * (ny+1) + 0;
}

struct Point
{
    Vector3 x;
    //double x;
    //double y;
    //double z;

    Point()
    {
        this->x(0) = 0.0;
        this->x(1) = 0.0;
        this->x(2) = 0.0;
    }

    Point(double x, double y, double z)
    {
        this->x(0) = x;
        this->x(1) = y;
        this->x(2) = z;
    }

    Point sum(double x, double y, double z)
    {
        Point c;

        c.x(0) = this->x(0) + x;
        c.x(1) = this->x(1) + y;
        c.x(2) = this->x(2) + z;

        return c;
    }
};

struct Line
{
    Point p0;
    Point p1;

    Line(Point p0, Point p1)
    {
        this->p0 = p0;
        this->p1 = p1;
    }

    double len()
    {
        Vector3 c = p1.x - p0.x;
        return c.len();
    }

    double dl(int n)
    {
        return len() / double(n);
    }

    double dx(int n)
    {
        return std::abs(p1.x(0) - p0.x(0)) / double(n);
    }

    double dy(int n)
    {
        return std::abs(p1.x(1) - p0.x(1)) / double(n);
    }

    double dz(int n)
    {
        return std::abs(p1.x(2) - p0.x(2)) / double(n);
    }

    std::vector<Point> generate(int n)
    {
        std::vector<Point> pts;

        pts.push_back(p0);

        for (int i=0; i<n-1; ++i)
        {
            pts.push_back(pts.back().sum(dx(n), dy(n), dz(n)));
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

struct Cell
{
    std::vector<IJK> ijk;

    Cell(std::vector<IJK> ijk)
    {
        this->ijk = ijk;
    }

    double volume()
    {
    }
};

std::vector<Cell> cells;

double triangle_area(Point& p0, Point& p1, Point& p2)
{
    return 0.5 * len(cross((p1.x - p0.x), (p2.x - p0.x)));
}

Vector3 vertex_centroid(std::vector<Point> points)
{
    Vector3 vc(0., 0., 0.);

    for (Point& p: points)
    {
        vc += p.x;
    }

    vc /= points.size();

    return vc;
}

void quad_area(std::vector<Point> points)
{
    area_ = 0.;
    Centroid vc = vertex_centroid(points_);
    PointPointer pvc = new Point(-1, vc);

    for (int i=0; i<points_.size(); ++i)
    {
        PointPointers pts;
        pts.push_back(points_[i]);

        if (i == points_.size() - 1)
        {
            pts.push_back(points_[0]);
        }
        else
        {
            pts.push_back(points_[i+1]);
        }

        pts.push_back(pvc);

        Area tri_area = triangle_area(pts);

        area_ += tri_area;
    }

    assert(area_ != 0.);
    assert(!std::isnan(area_));
}

void make_cells()
{
    for (int k=0; k<nz; ++k)
    {
        for (int j=0; j<ny; ++j)
        {
            for (int i=0; i<nx; ++i)
            {
                cells.push_back(IJK(i,j,k), IJK(i+1, j, k), IJK(i+1, j+1, k), IJK(i, j+1, k));
                cells.push_back(IJK(i,j,k+1), IJK(i+1, j, k+1), IJK(i+1, j+1, k+1), IJK(i, j+1, k+1));
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
            surfaces.push_back(IJK(i,j,0),  IJK(i+1, j, 0),  IJK(i+1, j+1, 0),  IJK(i, j+1, 0));
            surfaces.push_back(IJK(i,j,nz), IJK(i+1, j, nz), IJK(i+1, j+1, nz), IJK(i, j+1, nz));
        }
    }

    for (int j=0; j<ny; ++j)
    {
        for (int k=0; k<nz; ++k)
        {
            surfaces.push_back(IJK(0,j,k),  IJK(0, j, k+1),  IJK(0, j+1, k+1),  IJK(0, j+1, k));
            surfaces.push_back(IJK(nx,j,k),  IJK(nx, j, k+1),  IJK(nx, j+1, k+1),  IJK(nx, j+1, k));
        }
    }

    for (int i=0; i<nx; ++i)
    {
        for (int k=0; k<nz; ++k)
        {
            surfaces.push_back(IJK(i,0,k),  IJK(i+1, 0, k),  IJK(i+1, 0, k+1),  IJK(i, 0, k+1));
            surfaces.push_back(IJK(i,ny,k),  IJK(i+1, ny, k),  IJK(i+1, ny, k+1),  IJK(i, ny, k+1));
        }
    }
}

int main()
{
    Point a(0.0, 0.0, 0.0);
    Point b(1.0, 0.0, 0.0);
    Point c(1.0, 1.0, 0.0);
    Point d(0.0, 1.0, 0.0);

    Point e(0.0, 0.0, 1.0);
    Point f(1.0, 0.0, 1.0);
    Point g(1.0, 1.0, 1.0);
    Point h(0.0, 1.0, 1.0);

    Line l0(a, b);
    Line l1(b, c);
    Line l2(d, c);
    Line l3(a, d);

    Line l4(e, f);
    Line l5(f, g);
    Line l6(h, g);
    Line l7(e, h);

    Line l8(a, e);
    Line l9(b, f);
    Line l10(c, g);
    Line l11(d, h);

    Point p[nx+1][ny+1][nz+1];

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

    int ncell = nx * ny * nz;
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
                out << p[i][j][k].x(0) << " " << p[i][j][k].x(1) << " " << p[i][j][k].x(2) << std::endl;
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
                out << p[i][j][k].x(0) << " " << p[i][j][k].x(1) << " " << p[i][j][k].x(2) << std::endl;
            }
        }
    }
    out << "CELLS ";
    out << cells.size();
    out << " ";
    out << cells.size() * 9 + surfaces.size() * 5 << std::endl;
    for (Surface& surface: surfaces)
    {
        out << 4 << " ";
        for (IJK ijk: surface.ijk)
        {
            out << pos(ijk[0]) << " ";
            out << pos(ijk[1]) << " ";
            out << pos(ijk[2]) << " ";
            out << pos(ijk[3]) << " ";
        }
        out << std::endl;
    }

    for (Cell& cell: cells)
    {
        out << 8 << " ";
        for (IJK ijk: surface.ijk)
        {
            out << pos(ijk[0]) << " ";
            out << pos(ijk[1]) << " ";
            out << pos(ijk[2]) << " ";
            out << pos(ijk[3]) << " ";
            out << pos(ijk[4]) << " ";
            out << pos(ijk[5]) << " ";
            out << pos(ijk[6]) << " ";
            out << pos(ijk[7]) << " ";
        }

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

    //out << "CELL_DATA ";
    //out << ncell_sur << std::endl;
    //out << "SCALARS Area float 1" << std::endl;
    //out << "LOOKUP_TABLE default" << std::endl;
    //for (int i=0; i<ncell_xy; ++i)
    //{
    //        out << pos(i, j, 0) << " ";
    //        out << pos(i+1, j, 0) << " ";
    //        out << pos(i+1, j+1, 0) << " ";
    //        out << pos(i, j+1, 0) << " ";
    //    out << 9 << std::endl;
    //}







    out.close();

    return 0;
}
