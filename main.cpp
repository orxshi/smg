#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>

struct Point
{
    double x;
    double y;
    double z;

    Point()
    {
        this->x = 0.0;
        this->y = 0.0;
        this->z = 0.0;
    }

    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point sum(double x, double y, double z)
    {
        Point c;

        c.x = this->x + x;
        c.y = this->y + y;
        c.z = this->z + z;

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
        return std::sqrt(std::pow(p1.x - p0.x, 2.) + std::pow(p1.y - p0.y, 2.) + std::pow(p1.z - p0.z, 2.));
    }

    double dl(int n)
    {
        return len() / double(n);
    }

    double dx(int n)
    {
        return std::abs(p1.x - p0.x) / double(n);
    }

    double dy(int n)
    {
        return std::abs(p1.y - p0.y) / double(n);
    }

    double dz(int n)
    {
        return std::abs(p1.z - p0.z) / double(n);
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

    int nx = 2; // number of cells in x-direction
    int ny = 2;
    int nz = 2;

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

    for (int i=1; i<nx; ++i)
    {
        for (int k=0; k<nz+1; ++k)
        {
            Line l(p[i][0][k], p[i][ny][k]);
            auto pts = l.generate(ny);
            for (int j=1; j<ny; ++j)
            {
                p[i][j][k] = pts[j];
            }
        }
    }

    for (int k=1; k<nz; ++k)
    {
        for (int i=1; i<nx; ++i)
        {
            Line l(p[i][0][k], p[i][ny][k]);
            auto pts = l.generate(nz);
            for (int j=1; j<ny; ++j)
            {
                p[i][j][k] = pts[j];
            }
        }
    }

    for (int k=0; k<nz+1; ++k)
    {
        for (int i=0; i<nx+1; ++i)
        {
            for (int j=0; j<ny+1; ++j)
            {
                std::cout << i << "," << j << "," << k << " " << p[i][j][k].x << ", " << p[i][j][k].y << ", " << p[i][j][k].z << std::endl;
            }
        }
    }

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
    for (int i=0; i<nx+1; ++i)
    {
        for (int j=0; j<ny+1; ++j)
        {
		for (int k=0; k<nz+1; ++k)
		{
            out << p[i][j][k].x << " " << p[i][j][k].y << " " << p[i][j][k].z << std::endl;
		}
        }
    }

    out.close();

    return 0;
}
