#include "filters_par.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <iostream>
#include <string>
#include <cstdlib>

// to check for wrong usage of program.
int main(int argc, const char** argv)
{
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " [radius<=1000] [infile.ppm] [outfile.ppm] [num_threads]\n";
        return 1;
    }

    // to convert the entered cmd arguments into values to use later.
    int radius = std::stoi(argv[1]);
    std::string infile  = argv[2];
    std::string outfile = argv[3];
    int num_threads = std::stoi(argv[4]);

    if (radius < 0 || radius > 1000) {
        std::cerr << "Error: radius must be in [0,1000]\n";
        return 1;
    }
    if (num_threads < 1) num_threads = 1;

    try {
        //reads input image into matric object.
        PPM::Reader reader;
        Matrix m = reader(infile);
        //calls the pthreads blur function.
        Matrix out = Filter::blur_par(m, radius, num_threads);

        PPM::Writer writer;
        writer(out, outfile);
    } catch (const std::exception& e) {
        std::cerr << "Failed: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
