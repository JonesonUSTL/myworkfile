#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem {

// 每个节点统一采用 6 个自由度：平动 + 转动。
constexpr int kDofPerNode = 6;  // UX UY UZ RX RY RZ

// 支持的单元类型。
enum class ElementType { Truss2, Beam2, Shell4, Solid8 };

// 分析步类型：线性静力 / 非线性静力。
enum class AnalysisType { LinearStatic, NonlinearStatic };

// 材料本构类型。
enum class MaterialLaw { LinearElastic, BilinearElastoPlastic, J2Plasticity };

// 线性求解后端：致密高斯、Eigen 稀疏、PETSc、并行 CG。
enum class LinearSolverBackend { DenseGaussian, EigenSparse, PetscKsp, ParallelCG };

// 节点坐标定义。
struct Node {
  int id{};
  std::array<double, 3> x{0.0, 0.0, 0.0};
};

// 单元历史变量（当前为积分点等效塑性应变与应力占位）。
struct ElementState {
  std::vector<double> eqPlasticStrain;
  std::vector<std::array<double, 6>> stress;
};

// 通用单元定义：连接关系、截面、厚度、材料。
struct Element {
  int id{};
  ElementType type{ElementType::Truss2};
  std::vector<int> conn;
  double area{1.0};
  double thickness{1.0};
  std::string material;
  ElementState state;
};

// 材料参数定义。
struct Material {
  std::string name;
  MaterialLaw law{MaterialLaw::LinearElastic};
  double young{210e9};
  double poisson{0.3};
  double yieldStress{250e6};
  double hardening{1.0e9};
  double density{0.0};
  double conductivity{0.0};
  double specificHeat{0.0};
  double expansion{0.0};
  std::vector<std::pair<double, double>> plasticTable;  // (plastic_strain, stress)
};

// 位移边界条件：支持节点号或节点集名称。
struct DofBC {
  int nodeId{};              // >0 时使用节点号
  std::string nodeSetName;   // 非空时使用节点集
  int dof{};                 // 1..6
  double value{};
};

// 节点载荷：支持节点号或节点集名称。
struct NodalLoad {
  int nodeId{};              // >0 时使用节点号
  std::string nodeSetName;   // 非空时使用节点集
  int dof{};                 // 1..6
  double value{};
  std::string amplitudeName;
};

// 体载荷/面载荷的轻量抽象（当前用于 DLOAD/DSLOAD 统一处理）。
struct BodyLoad {
  std::string elset;
  std::string type;          // 如 P / PX / PY / PZ / GRAV
  double value{};
  std::string amplitudeName;
};

// 幅值曲线（时间-系数）。
struct Amplitude {
  std::string name;
  std::vector<std::pair<double, double>> points;
};

// 耦合约束（当前保存解析结果，后续可扩展到方程约束组装）。
struct Coupling {
  int referenceNode{};
  std::string surfaceName;
  std::vector<int> surfaceNodes;
  std::vector<int> dofs;
  double penalty{1e10};
};

// 接触对（法向 + 摩擦罚函数参数）。
struct ContactPair {
  std::string masterSurface;
  std::string slaveSurface;
  double penalty{1e9};
  double friction{0.3};
};


// 多点约束方程（MPC）：slaveDof = sum(coeff_i * masterDof_i) + offset。
struct MpcConstraint {
  int slaveNode{};
  int slaveDof{};
  std::vector<int> masterNodes;
  std::vector<int> masterDofs;
  std::vector<double> coefficients;
  double offset{0.0};
};

// 输出请求（当前保存关键字，便于后续做选择性输出）。
struct OutputRequest {
  std::string kind;  // NODE OUTPUT / ELEMENT OUTPUT
  std::vector<std::string> variables;
};

// 分析步参数。
struct Step {
  std::string name{"Step-1"};
  AnalysisType type{AnalysisType::LinearStatic};
  int increments{10};
  double totalLoadFactor{1.0};
  double initialIncrement{0.1};
  double minIncrement{1e-6};
  double maxIncrement{0.25};
  int maxNewtonIters{30};
  double tolerance{1e-8};
  int maxLinearIters{5000};
  int maxCutbacks{8};
  LinearSolverBackend solver{LinearSolverBackend::DenseGaussian};
  std::string amplitudeName;
  bool useArcLength{false};
  double arcLengthRadius{1e-2};
  double arcLengthMinRadius{1e-5};
  double arcLengthMaxRadius{5e-2};
  double arcLengthGrowFactor{1.25};
  double arcLengthShrinkFactor{0.5};
  int fieldOutputFrequency{1};
};



// 热力耦合扩展接口（当前为预留数据结构，便于后续热-结构耦合实现）。
struct TemperatureBC {
  int nodeId{};
  std::string nodeSetName;
  double temperature{};
};

struct HeatLoad {
  int nodeId{};
  std::string nodeSetName;
  double value{};
};

// 有限元模型容器。
struct Model {
  std::string modelName{"Model-1"};
  std::vector<std::string> partNames;
  std::vector<Node> nodes;
  std::vector<Element> elements;
  std::unordered_map<std::string, Material> materials;
  std::unordered_map<std::string, std::vector<int>> nsets;
  std::unordered_map<std::string, std::vector<int>> elsets;
  std::unordered_map<std::string, Amplitude> amplitudes;
  std::vector<Coupling> couplings;
  std::vector<ContactPair> contacts;
  std::vector<MpcConstraint> mpcs;
  std::vector<OutputRequest> outputRequests;
  std::vector<DofBC> bcs;
  std::vector<NodalLoad> loads;
  std::vector<BodyLoad> bodyLoads;
  std::vector<TemperatureBC> temperatureBCs;
  std::vector<HeatLoad> heatLoads;
  Step step;
};


// 计算结果容器。
struct Result {
  std::vector<double> displacement;
  std::vector<double> reaction;
  std::vector<double> elementAxialForce;
  std::vector<double> elementStrain;
  std::vector<double> elementStress;
  std::vector<double> elementVonMises;
  std::vector<std::vector<double>> displacementFrames;
  std::vector<std::vector<double>> reactionFrames;
  std::vector<double> frameLoadFactors;
};

}  // namespace fem
