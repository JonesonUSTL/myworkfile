#include "vtk_writer.hpp"

#include <fstream>
#include <stdexcept>

namespace fem {
namespace {

int nodeIndex(const Model& model, int nodeId) {
  for (size_t i = 0; i < model.nodes.size(); ++i) {
    if (model.nodes[i].id == nodeId) return static_cast<int>(i);
  }
  throw std::runtime_error("Element references unknown node");
}

int vtkCellType(ElementType t) {
  if (t == ElementType::Truss2 || t == ElementType::Beam2) return 3;   // line
  if (t == ElementType::Shell4) return 9;                                // quad
  if (t == ElementType::Solid8) return 12;                               // hexahedron
  return 3;
}

}  // namespace

void writeVTK(const std::string& filePath, const Model& model, const Result& result, double scale) {
  std::ofstream out(filePath);
  if (!out) throw std::runtime_error("Cannot write VTK file: " + filePath);

  out << "# vtk DataFile Version 3.0\n";
  out << "FEM result\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  out << "POINTS " << model.nodes.size() << " double\n";
  for (size_t i = 0; i < model.nodes.size(); ++i) {
    const auto& n = model.nodes[i];
    out << n.x[0] + result.displacement[i * kDofPerNode + 0] * scale << " "
        << n.x[1] + result.displacement[i * kDofPerNode + 1] * scale << " "
        << n.x[2] + result.displacement[i * kDofPerNode + 2] * scale << "\n";
  }

  size_t cellSize = 0;
  for (const auto& e : model.elements) cellSize += 1 + e.conn.size();
  out << "CELLS " << model.elements.size() << " " << cellSize << "\n";
  for (const auto& e : model.elements) {
    out << e.conn.size();
    for (int nid : e.conn) out << " " << nodeIndex(model, nid);
    out << "\n";
  }

  out << "CELL_TYPES " << model.elements.size() << "\n";
  for (const auto& e : model.elements) out << vtkCellType(e.type) << "\n";

  out << "POINT_DATA " << model.nodes.size() << "\n";
  out << "VECTORS displacement double\n";
  for (size_t i = 0; i < model.nodes.size(); ++i) {
    out << result.displacement[i * kDofPerNode + 0] << " " << result.displacement[i * kDofPerNode + 1] << " "
        << result.displacement[i * kDofPerNode + 2] << "\n";
  }

  out << "VECTORS reaction double\n";
  for (size_t i = 0; i < model.nodes.size(); ++i) {
    out << result.reaction[i * kDofPerNode + 0] << " " << result.reaction[i * kDofPerNode + 1] << " "
        << result.reaction[i * kDofPerNode + 2] << "\n";
  }

  out << "CELL_DATA " << model.elements.size() << "\n";
  out << "SCALARS axial_force double 1\nLOOKUP_TABLE default\n";
  for (double v : result.elementAxialForce) out << v << "\n";

  out << "SCALARS strain double 1\nLOOKUP_TABLE default\n";
  for (double v : result.elementStrain) out << v << "\n";

  out << "SCALARS stress double 1\nLOOKUP_TABLE default\n";
  for (double v : result.elementStress) out << v << "\n";

  out << "SCALARS von_mises double 1\nLOOKUP_TABLE default\n";
  for (double v : result.elementVonMises) out << v << "\n";
}

}  // namespace fem
