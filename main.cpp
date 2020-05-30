#include <iostream>
#include <mpi.h>
#include <boost/program_options.hpp>
#include <fstream>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    int commsize, rank, numprocesses;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&commsize);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    numprocesses = commsize-1;

    if (rank == 0) {
        std::string conf_file;
        if (argc == 1) {
            conf_file = "config.dat";
        } else if (argc == 2) {
            conf_file = argv[1];
        }

        std::ifstream conf(conf_file, std::ios::in);

        if (!conf.is_open()) {
            std::cerr << "Could not open the configuration file. Set your working directory to ..\n";
            exit(2);
        }

        double cap, conduct, dens;
        int x, y, delta_x, delta_y, delta_t, iter, save_step;

        po::options_description config_parser;
        config_parser.add_options()
                ("capacity", po::value<double>(&cap))
                ("conductivity", po::value<double>(&conduct))
                ("density", po::value<double>(&dens))
                ("x", po::value<int>(&x))
                ("y", po::value<int>(&y))
                ("delta_x", po::value<int>(&delta_x)->default_value(1))
                ("delta_y", po::value<int>(&delta_y)->default_value(1))
                ("delta_t", po::value<int>(&delta_t))
                ("iterations", po::value<int>(&iter))
                ("save_step", po::value<int>(&save_step));
        po::variables_map vm;
        store(parse_config_file(conf, config_parser), vm);
        notify(vm);

        std::cout << cap << std::endl;
        MPI_Finalize();
    }

    if (rank != 0) {
        std::cout << rank << std::endl;
        MPI_Finalize();
    }



    return 0;
}
