cmake_minimum_required(VERSION 3.16)
project(mini_fem_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FEM_ENABLE_EIGEN "Enable Eigen sparse linear solver backend" ON)
option(FEM_ENABLE_PETSC "Enable PETSc KSP backend" ON)
option(FEM_ENABLE_OPENMP "Enable OpenMP for parallel assembly/solver" ON)

add_executable(fem_solver
    src/main.cpp
    src/inp_parser.cpp
    src/solver.cpp
    src/vtk_writer.cpp
)

target_include_directories(fem_solver PRIVATE include)

if (FEM_ENABLE_EIGEN)
  find_package(Eigen3 QUIET)
  if (Eigen3_FOUND)
    target_link_libraries(fem_solver PRIVATE Eigen3::Eigen)
    target_compile_definitions(fem_solver PRIVATE FEM_USE_EIGEN)
    message(STATUS "Eigen backend enabled")
  else()
    message(STATUS "Eigen not found, fallback to dense Gaussian solver")
  endif()
endif()

if (FEM_ENABLE_PETSC)
  find_package(PETSc QUIET)
  if (PETSc_FOUND)
    target_include_directories(fem_solver PRIVATE ${PETSC_INCLUDES})
    target_link_libraries(fem_solver PRIVATE ${PETSC_LIBRARIES})
    target_compile_definitions(fem_solver PRIVATE FEM_USE_PETSC)
    message(STATUS "PETSc backend enabled")
  else()
    message(STATUS "PETSc not found, petsc backend unavailable at runtime")
  endif()
endif()

if (FEM_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(fem_solver PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled")
  else()
    message(STATUS "OpenMP not found, run in serial mode")
  endif()
endif()

if (MSVC)
  target_compile_options(fem_solver PRIVATE /W4)
else()
  target_compile_options(fem_solver PRIVATE -Wall -Wextra -Wpedantic)
endif()
