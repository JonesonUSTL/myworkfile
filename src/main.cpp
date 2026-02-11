#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "inp_parser.hpp"
#include "solver.hpp"
#include "vtk_writer.hpp"

int main(int argc, char** argv) {
  try {
    if (argc < 3) {
      std::cerr << "Usage: fem_solver <input.inp> <output.vtk> [deform_scale] [solver=dense|eigen|petsc|pcg]\n";
      return 1;
    }

    const std::string inp = argv[1];
    const std::string vtk = argv[2];
    const double scale = (argc >= 4) ? std::stod(argv[3]) : 1.0;

    auto model = fem::parseInpFile(inp);
    if (argc >= 5) {
      std::string s = argv[4];
      if (s == "dense") model.step.solver = fem::LinearSolverBackend::DenseGaussian;
      else if (s == "eigen") model.step.solver = fem::LinearSolverBackend::EigenSparse;
      else if (s == "petsc") model.step.solver = fem::LinearSolverBackend::PetscKsp;
      else if (s == "pcg") model.step.solver = fem::LinearSolverBackend::ParallelCG;
    }

    auto result = fem::solve(model);

    const auto outPath = std::filesystem::path(vtk);
    if (outPath.has_parent_path()) std::filesystem::create_directories(outPath.parent_path());
    fem::writeVTK(vtk, model, result, scale);

    std::cout << "Solve completed.\n";
    std::cout << "Nodes: " << model.nodes.size() << ", Elements: " << model.elements.size() << "\n";
    std::cout << "Analysis: "
              << (model.step.type == fem::AnalysisType::NonlinearStatic ? "Nonlinear static" : "Linear static")
              << "\n";
    std::cout << "VTK output: " << vtk << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}
