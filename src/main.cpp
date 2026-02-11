#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "inp_parser.hpp"
#include "solver.hpp"
#include "vtk_writer.hpp"

namespace {

// 命令行帮助文本：覆盖安装、运行、示例与求解器后端说明。
void printHelp(const char* exe) {
  std::cout
      << "\n==== mini_fem_solver 使用说明 ====\n"
      << "用法:\n"
      << "  " << exe
      << " <input.inp> <output.vtk> [deform_scale] [solver=dense|eigen|petsc|pcg]\n\n"
      << "常用命令:\n"
      << "  " << exe << " --help                  # 查看帮助\n"
      << "  " << exe << " --examples              # 列出内置示例\n"
      << "  " << exe << " examples/linear_truss.inp output/linear.vtk 1.0 dense\n\n"
      << "参数说明:\n"
      << "  input.inp      Abaqus 风格输入文件\n"
      << "  output.vtk     结果输出路径（会自动创建目录）\n"
      << "  deform_scale   变形显示缩放，默认 1.0\n"
      << "  solver         dense/eigen/petsc/pcg，默认 dense\n\n"
      << "建议阅读文档:\n"
      << "  README.md\n"
      << "  docs/user_guide_cn.md\n"
      << "===============================\n\n";
}

// 打印仓库中推荐入门案例，便于快速试算。
void printExamples() {
  const std::vector<std::string> demos = {
      "examples/linear_truss.inp       # 线性桁架（最小入门）",
      "examples/nonlinear_truss.inp    # 非线性桁架（Newton-Raphson）",
      "examples/truss_plastic.inp      # 双线性弹塑性桁架",
      "examples/c3d8_linear.inp        # C3D8 实体线性静力",
      "examples/classic_bar_tension.inp# 经典单杆受拉",
      "examples/classic_beam_tip_load.inp # 经典悬臂梁端载",
      "examples/classic_solid_patch.inp# 经典实体 patch/压缩算例",
      "examples/advanced_keywords.inp  # 接触/耦合/MPC/RIKS 扩展关键字",
  };
  std::cout << "\n内置示例列表:\n";
  for (const auto& e : demos) std::cout << "  - " << e << "\n";
  std::cout << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    // 交互入口：无参数或帮助参数时直接给出说明，提升可用性。
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
      printHelp(argv[0]);
      return (argc < 2) ? 1 : 0;
    }
    if (std::string(argv[1]) == "--examples") {
      printExamples();
      return 0;
    }

    if (argc < 3) {
      std::cerr << "参数不足，请使用 --help 查看完整说明。\n";
      return 1;
    }

    const std::string inp = argv[1];
    const std::string vtk = argv[2];
    const double scale = (argc >= 4) ? std::stod(argv[3]) : 1.0;

    // 1) 解析输入模型。
    auto model = fem::parseInpFile(inp);

    // 2) 解析用户指定的线性求解器后端。
    if (argc >= 5) {
      std::string s = argv[4];
      if (s == "dense") model.step.solver = fem::LinearSolverBackend::DenseGaussian;
      else if (s == "eigen") model.step.solver = fem::LinearSolverBackend::EigenSparse;
      else if (s == "petsc") model.step.solver = fem::LinearSolverBackend::PetscKsp;
      else if (s == "pcg") model.step.solver = fem::LinearSolverBackend::ParallelCG;
      else {
        std::cerr << "未知求解器后端: " << s << "，将回退到 dense。\n";
      }
    }

    // 3) 执行求解。
    auto result = fem::solve(model);

    // 4) 输出 VTK 结果。
    const auto outPath = std::filesystem::path(vtk);
    if (outPath.has_parent_path()) std::filesystem::create_directories(outPath.parent_path());
    fem::writeVTK(vtk, model, result, scale);

    // 5) 控制台摘要，便于批处理日志回看。
    std::cout << "Solve completed.\n";
    std::cout << "Nodes: " << model.nodes.size() << ", Elements: " << model.elements.size() << "\n";
    std::cout << "Analysis: "
              << (model.step.type == fem::AnalysisType::NonlinearStatic ? "Nonlinear static" : "Linear static")
              << "\n";
    std::cout << "Displacement size: " << result.displacement.size() << ", Reaction size: " << result.reaction.size()
              << "\n";
    std::cout << "VTK output: " << vtk << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 2;
  }
}
