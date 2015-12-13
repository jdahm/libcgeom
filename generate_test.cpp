#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <stdexcept>

constexpr double pi = 3.141592653589793238463;

void generate_regular(const std::string& fname, unsigned int narg, char *arg[])
{
        if (narg != 7) throw std::runtime_error("USAGE: xstart dx nx ystart dy ny n");

        double xs, dx, ys, dy;
        unsigned int nx, ny, n;
        { std::stringstream ss(arg[0]); ss >> xs; }
        { std::stringstream ss(arg[1]); ss >> dx; }
        { std::stringstream ss(arg[2]); ss >> nx; }
        { std::stringstream ss(arg[3]); ss >> ys; }
        { std::stringstream ss(arg[4]); ss >> dy; }
        { std::stringstream ss(arg[5]); ss >> ny; }
        { std::stringstream ss(arg[6]); ss >> n; }

        if (nx < 2 || ny < 2) throw std::runtime_error("More perimeter points needed");

        std::ofstream fst(fname, std::ios::out);

        const unsigned int np = nx * ny;
        fst << np + n << " " << "0" << std::endl;

        for (unsigned int i=0; i<nx; i++)
                for (unsigned int j=0; j<ny; j++)
                        fst << xs + i*dx << " " << ys + j*dy << std::endl;

        const double xe = xs + (nx - 1) * dx;
        const double ye = ys + (ny - 1) * dy;

        std::random_device rd;
        std::mt19937 gen(rd());

        for (unsigned int i=0; i<n; i++) {
                const double x = xs + (xe - xs) * std::generate_canonical<double, 512>(gen);
                const double y = ys + (ye - ys) * std::generate_canonical<double, 512>(gen);
                fst << x << " " << y << std::endl;
        }

        fst.close();
}

unsigned int perimeter_density(unsigned int np)
{
        // Uses heuristics to calculate density of points on perimeter to "look good"
        return std::min(
                0.2 * np < 2 ? 2 : 0.2 * np,
                std::max(2.0,
                         0.25 * std::sqrt(static_cast<double>(np)))
                );
}

std::vector<double> perimeter(double xs, double xe, double ys, double ye,
                            unsigned int nx, unsigned int ny)
{
        std::vector<double> point;

        if (nx < 2 || ny < 2) throw std::runtime_error("More perimeter points needed");

        const double dy = (ye - ys) / (ny - 1);
        const double dx = (xe - xs) / (nx - 1);

        // Left
        {
                double x = xs;
                double y = ys;
                while (y <= ye) {
                        point.push_back(x);
                        point.push_back(y);
                        y += dy;
                }
        }
        // Top
        {
                double x = xs + dx;
                double y = ye;
                while (x <= xe) {
                        point.push_back(x);
                        point.push_back(y);
                        x += dx;
                }
        }
        // Right
        {
                double x = xe;
                double y = ye - dy;
                while (y >= ys) {
                        point.push_back(x);
                        point.push_back(y);
                        y -= dy;
                }
        }
        // Bottom
        {
                double x = xe - dx;
                double y = ys;
                while (x > xs) {
                        point.push_back(x);
                        point.push_back(y);
                        x -= dx;
                }
        }

        return point;
}

void generate_uniform(const std::string& fname, unsigned int narg, char *arg[])
{
        if (narg != 6) throw std::runtime_error("USAGE: xstart xend ystart yend np use_perim");

        double xs, xe, ys, ye;
        unsigned int np;
        bool use_perim;

        { std::stringstream ss(arg[0]); ss >> xs; }
        { std::stringstream ss(arg[1]); ss >> xe; }
        { std::stringstream ss(arg[2]); ss >> ys; }
        { std::stringstream ss(arg[3]); ss >> ye; }
        { std::stringstream ss(arg[4]); ss >> np; }
        { std::stringstream ss(arg[5]); ss >> use_perim; }

        unsigned int nperim = 0;
        std::vector<double> perim;
        if (use_perim) {
                const unsigned int per = perimeter_density(np);
                perim = perimeter(0, 1, 0, 1, per, per);
                nperim = perim.size() / 2;
        }

        if (nperim > np) throw std::runtime_error("nperim > np");

        std::ofstream fst(fname, std::ios::out);

        fst << np << " " << "0" << std::endl;

        if (use_perim)
                for (unsigned int i=0; i<perim.size(); i+=2)
                        fst << perim[i] << " " << perim[i+1] << std::endl;

        std::random_device rd;
        std::mt19937 gen(rd());

        for (unsigned int i=0; i<np - nperim; i++) {
                const double x = xs + (xe - xs) * std::generate_canonical<double, 512>(gen);
                const double y = ys + (ye - ys) * std::generate_canonical<double, 512>(gen);
                fst << x << " " << y << std::endl;
        }

        fst.close();
}

void generate_gaussian(const std::string& fname, unsigned int narg, char *arg[])
{
        if (narg != 9) throw std::runtime_error("USAGE: xstart xend ystart yend xc yc sigma np use_perim");

        double xs, xe, ys, ye, xc, yc, sigma;
        unsigned int np;
        bool use_perim;

        { std::stringstream ss(arg[0]); ss >> xs; }
        { std::stringstream ss(arg[1]); ss >> xe; }
        { std::stringstream ss(arg[2]); ss >> ys; }
        { std::stringstream ss(arg[3]); ss >> ye; }
        { std::stringstream ss(arg[4]); ss >> xc; }
        { std::stringstream ss(arg[5]); ss >> yc; }
        { std::stringstream ss(arg[6]); ss >> sigma; }
        { std::stringstream ss(arg[7]); ss >> np; }
        { std::stringstream ss(arg[8]); ss >> use_perim; }

        unsigned int nperim = 0;
        std::vector<double> perim;
        if (use_perim) {
                const unsigned int per = perimeter_density(np);
                perim = perimeter(0, 1, 0, 1, per, per);
                nperim = perim.size() / 2;
        }

        if (nperim > np) throw std::runtime_error("nperim > np");

        std::ofstream fst(fname, std::ios::out);
        fst << np << " " << "0" << std::endl;

        if (use_perim)
                for (unsigned int i=0; i<perim.size(); i+=2)
                        fst << perim[i] << " " << perim[i+1] << std::endl;

        std::random_device rd;
        std::mt19937 gen(rd());

        std::vector<double> b{0, sigma, 1};
        std::vector<double> w{1, 0.05, 0};
        std::piecewise_linear_distribution<double> d(b.begin(), b.end(), w.begin());
        // std::normal_distribution<> d(0, sigma);
        std::cout << "np = " << np << std::endl;
        unsigned int i=0;
        while (i < np - nperim) {
                const double dist = 2 * d(gen);
                const double angle = 2 * pi * std::generate_canonical<double, 512>(gen);
                const double x = xc + dist * std::cos(angle);
                const double y = yc + dist * std::sin(angle);
                std::cout << x << " " << y << std::endl;
                if ((x > xs) && (x < xe) && (y > ys) && (y < ye)) {
                        fst << x << " " << y << std::endl;
                        i++;
                }
        }

        fst.close();
}

void generate_line(const std::string& fname, unsigned int narg, char *arg[])
{
        if (narg != 3) throw std::runtime_error("USAGE: sigma np use_perim");

        double sigma;
        unsigned int np;
        bool use_perim;

        { std::stringstream ss(arg[0]); ss >> sigma; }
        { std::stringstream ss(arg[1]); ss >> np; }
        { std::stringstream ss(arg[2]); ss >> use_perim; }

        unsigned int nperim = 0;
        std::vector<double> perim;
        if (use_perim) {
                const unsigned int per = perimeter_density(np);
                perim = perimeter(0, 1, 0, 1, per, per);
                nperim = perim.size() / 2;
        }

        if (nperim > np) throw std::runtime_error("nperim > np");

        std::ofstream fst(fname, std::ios::out);
        fst << np << " " << "0" << std::endl;

        std::random_device rd;
        std::mt19937 gen(rd());

        if (use_perim)
                for (unsigned int i=0; i<perim.size(); i+=2)
                        fst << perim[i] << " " << perim[i+1] << std::endl;

        std::vector<double> b{0, sigma, 1};
        std::vector<double> w{1, 0.05, 0};
        std::piecewise_linear_distribution<double> d(b.begin(), b.end(), w.begin());
        // std::normal_distribution<> d(0.5, sigma);
        unsigned int i=0;
        while (i < np - nperim) {
                const double x = (i%2) ? 0.5*(1 + d(gen)) : 0.5*(1 - d(gen));
                const double y = std::generate_canonical<double, 512>(gen);
                fst << x << " " << y << std::endl;
                i++;
        }

        fst.close();
}

void generate_circle(const std::string& fname, unsigned int narg, char *arg[])
{
        if (narg != 1) throw std::runtime_error("USAGE: np");

        unsigned int np;

        { std::stringstream ss(arg[0]); ss >> np; }

        std::ofstream fst(fname, std::ios::out);
        fst << np << " " << "0" << std::endl;

        const double dalpha = 2 * pi / np;
        for (unsigned int i=0; i<np; i++) {
                const double x = cos(i * dalpha);
                const double y = sin(i * dalpha);
                fst << x << " " << y << std::endl;
        }

        fst.close();
}

int main(int argc, char *argv[]) {
        if (argc < 3)
                throw std::runtime_error("USAGE: ./generate_test filename distribution (other options)");

        const std::string fname(argv[1]);
        const std::string s(argv[2]);

        const unsigned int offset = 3;
        const unsigned int narg = argc - offset;
        if (s == "regular") generate_regular(fname, narg, argv + offset);
        else if (s == "uniform") generate_uniform(fname, narg, argv + offset);
        else if (s == "gaussian") generate_gaussian(fname, narg, argv + offset);
        else if (s == "line") generate_line(fname, narg, argv + offset);
        else if (s == "circle") generate_circle(fname, narg, argv + offset);
        else
                throw std::runtime_error("Unknown distribution");

        return 0;
}

