cmake_minimum_required(VERSION 3.15)
project(mpi_heat)

set(CMAKE_CXX_STANDARD 20)

find_package(MPI REQUIRED)

find_package(Boost REQUIRED COMPONENTS
        locale filesystem system program_options iostreams)
find_package(ImageMagick COMPONENTS Magick++)
include_directories(${ImageMagick_INCLUDE_DIRS})
set(Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
set(Boost_USE_MULTITHREAD ON)  # enable multithreading

include_directories(${Boost_INCLUDE_DIRS})
add_executable(mpi_heat main.cpp)
target_include_directories(mpi_heat PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_compile_options(mpi_heat PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(mpi_heat ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${Boost_LIBRARIES} ${ImageMagick_LIBRARIES})
