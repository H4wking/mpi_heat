#include <iostream>
#include <mpi.h>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
#include <fstream>
#include <mutex>

namespace po = boost::program_options;

std::mutex m;

int main(int argc, char *argv[]) {
    int commsize, rank, numprocesses;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&commsize);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    numprocesses = commsize-1;

    if (rank == 0) {
        std::string conf_file;
        if (argc == 1) {
            conf_file = "../config.dat";
        } else if (argc == 2) {
            conf_file = argv[1];
        }

        std::ifstream conf(conf_file, std::ios::in);

        if (!conf.is_open()) {
            std::cerr << "Could not open the configuration file. Set your working directory to ..\n";
            exit(2);
        }

        double cap, conduct, dens, delta_t;
        int x, y, delta_x, delta_y, iter, save_step;

        po::options_description config_parser;
        config_parser.add_options()
                ("heat_capacity", po::value<double>(&cap))
                ("thermal_conductivity", po::value<double>(&conduct))
                ("density", po::value<double>(&dens))
                ("x", po::value<int>(&x))
                ("y", po::value<int>(&y))
                ("delta_x", po::value<int>(&delta_x)->default_value(1))
                ("delta_y", po::value<int>(&delta_y)->default_value(1))
                ("delta_t", po::value<double>(&delta_t))
                ("iterations", po::value<int>(&iter))
                ("save_step", po::value<int>(&save_step));
        po::variables_map vm;
        store(parse_config_file(conf, config_parser), vm);
        notify(vm);

        typedef boost::multi_array<double, 2> array_type;
        array_type A(boost::extents[y][x]);

        std::ifstream matrix_file {"../matrix.txt"};
        if (!matrix_file.is_open()) return -1;

        for (int i = 0; i < y; i++) {
            for (int j = 0; j < x; j++) {
                matrix_file >> A[i][j];
            }
        }

//        for (int i = 0; i < y; i++) {
//            for (int j = 0; j < x; j++) {
//                std::cout << A[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }

        double alpha = conduct / (dens * cap);

        int avrows = y / numprocesses;
        int extra = y % numprocesses;

        int rows;
        int offset = 0;

        for (int i = 1; i <= numprocesses; i++) {
            rows = (i <= extra) ? avrows + 1 : avrows;

//            if (i == 1 || i == numprocesses) {
//                rows++;
//            } else {
//                rows += 2;
//            }

            MPI_Send(&alpha, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(&x, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&delta_t, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(&iter, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&save_step, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&A[offset][0], rows * x, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

            offset += rows;
        }

//        std::cout << alpha << std::endl;
        MPI_Finalize();
    }

    if (rank != 0) {
        double alpha, delta_t;
        int x, rows, iter, save_step;
        MPI_Status status;

        MPI_Recv(&alpha, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&delta_t, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&iter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&save_step, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        typedef boost::multi_array<double, 2> array_type;

        if (rank == 1 || rank == numprocesses) {
            rows++;
        } else {
            rows += 2;
        }

        array_type Ai(boost::extents[rows][x]);

        if (rank == 1) {
            MPI_Recv(&Ai[0][0], rows * x, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&Ai[1][0], rows * x, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }

//        usleep(10000 * rank);
//
//        std::cout << rank << std::endl;
//
//        for (int i = 0; i < rows; i++) {
//            for (int j = 0; j < x; j++) {
//                std::cout << Ai[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }
//
//
//        std::cout << std::endl;

        for (int i = 0; i < 1; i++) {
            if (rank != numprocesses) {
                MPI_Sendrecv(&Ai[rows-2][0], x, MPI_DOUBLE, rank + 1, 1,
                             &Ai[rows-1][0], x, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            }
            if (rank != 1) {
                MPI_Sendrecv(&Ai[1][0], x, MPI_DOUBLE, rank - 1, 1,
                             &Ai[0][0], x, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            }
        }

        usleep(10000 * rank);

        std::cout << rank << std::endl;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < x; j++) {
                std::cout << Ai[i][j] << " ";
            }
            std::cout << std::endl;
        }


        std::cout << std::endl;


        MPI_Finalize();
    }

    return 0;
}
