#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef FEM_USE_EIGEN
#include <Eigen/Sparse>
#endif
#ifdef FEM_USE_PETSC
#include <petscksp.h>
#endif

namespace fem {
namespace {

using Matrix = std::vector<std::vector<double>>;

int nodeIdx(const Model& model, int nodeId) {
  for (size_t i = 0; i < model.nodes.size(); ++i)
    if (model.nodes[i].id == nodeId) return static_cast<int>(i);
  throw std::runtime_error("Unknown node id");
}


// 根据节点号或节点集展开目标节点列表。
std::vector<int> expandNodes(const Model& model, int nodeId, const std::string& setName) {
  std::vector<int> ids;
  if (nodeId > 0) {
    ids.push_back(nodeId);
  } else if (!setName.empty()) {
    auto it = model.nsets.find(setName);
    if (it != model.nsets.end()) ids = it->second;
  }
  return ids;
}

// 计算向量点积（支持 OpenMP 并行归约）。
double dotVec(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : s)
#endif
  for (int i = 0; i < static_cast<int>(a.size()); ++i) s += a[i] * b[i];
  return s;
}

// 矩阵向量乘（致密存储下采用 OpenMP 行并行）。
std::vector<double> matVec(const Matrix& A, const std::vector<double>& x) {
  std::vector<double> y(x.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < static_cast<int>(A.size()); ++i) {
    double sum = 0.0;
    for (int j = 0; j < static_cast<int>(x.size()); ++j) sum += A[i][j] * x[j];
    y[i] = sum;
  }
  return y;
}

// 并行 CG 迭代求解（SPD 假设，配合正则化和边界后方程可用）。
std::vector<double> parallelCgSolve(const Matrix& A, const std::vector<double>& b, int maxIter, double tol) {
  const int n = static_cast<int>(b.size());
  std::vector<double> x(n, 0.0), r = b, p = r;
  double rsold = dotVec(r, r);
  if (std::sqrt(rsold) < tol) return x;
  for (int it = 0; it < maxIter; ++it) {
    auto Ap = matVec(A, p);
    double denom = dotVec(p, Ap);
    if (std::abs(denom) < 1e-30) break;
    double alpha = rsold / denom;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    }
    double rsnew = dotVec(r, r);
    if (std::sqrt(rsnew) < tol) break;
    double beta = rsnew / rsold;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];
    rsold = rsnew;
  }
  return x;
}


// 前向声明：耦合/接触罚函数。
void applyCouplingPenalty(const Model& model, Matrix& K, std::vector<double>* fint, const std::vector<double>* u);
void applyContactPenalty(const Model& model, Matrix& K, std::vector<double>* fint, const std::vector<double>* u);
std::vector<double> solveWithMpcLagrange(const Model& model, const Matrix& K, const std::vector<double>& rhs, const std::vector<double>* uCurrent);

std::array<double, 3> unitVec(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (L < 1e-20) throw std::runtime_error("Zero-length element");
  return {dx / L, dy / L, dz / L};
}

std::vector<double> gaussianSolve(Matrix A, std::vector<double> b) {
  const int n = static_cast<int>(b.size());
  for (int i = 0; i < n; ++i) {
    int p = i;
    for (int r = i + 1; r < n; ++r)
      if (std::abs(A[r][i]) > std::abs(A[p][i])) p = r;
    if (std::abs(A[p][i]) < 1e-20) throw std::runtime_error("Singular matrix");
    if (p != i) {
      std::swap(A[p], A[i]);
      std::swap(b[p], b[i]);
    }
    const double d = A[i][i];
    for (int c = i; c < n; ++c) A[i][c] /= d;
    b[i] /= d;
    for (int r = 0; r < n; ++r) {
      if (r == i) continue;
      const double f = A[r][i];
      if (std::abs(f) < 1e-30) continue;
      for (int c = i; c < n; ++c) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  return b;
}

std::vector<double> solveLinearSystem(const Model& model, const Matrix& K, const std::vector<double>& rhs) {
#ifdef FEM_USE_PETSC
  if (model.step.solver == LinearSolverBackend::PetscKsp) {
    static bool petscInit = false;
    if (!petscInit) {
      int argc = 0;
      char** argv = nullptr;
      PetscInitialize(&argc, &argv, nullptr, nullptr);
      petscInit = true;
    }

    const PetscInt n = static_cast<PetscInt>(rhs.size());
    Mat A;
    Vec b, x;
    KSP ksp;

    MatCreate(PETSC_COMM_SELF, &A);
    MatSetSizes(A, n, n, n, n);
    MatSetFromOptions(A);
    MatSetUp(A);
    for (PetscInt i = 0; i < n; ++i) {
      for (PetscInt j = 0; j < n; ++j) {
        const double v = K[i][j];
        if (std::abs(v) > 0.0) MatSetValue(A, i, j, v, INSERT_VALUES);
      }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    VecCreateSeq(PETSC_COMM_SELF, n, &b);
    VecDuplicate(b, &x);
    for (PetscInt i = 0; i < n; ++i) VecSetValue(b, i, rhs[i], INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPCG);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCGAMG);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, b, x);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    if (reason < 0) throw std::runtime_error("PETSc KSP diverged");

    std::vector<double> out(n, 0.0);
    const PetscScalar* arr = nullptr;
    VecGetArrayRead(x, &arr);
    for (PetscInt i = 0; i < n; ++i) out[i] = arr[i];
    VecRestoreArrayRead(x, &arr);

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    return out;
  }
#endif

#ifdef FEM_USE_EIGEN
  if (model.step.solver == LinearSolverBackend::EigenSparse) {
    using SpMat = Eigen::SparseMatrix<double>;
    std::vector<Eigen::Triplet<double>> tri;
    const int n = static_cast<int>(rhs.size());
    tri.reserve(static_cast<size_t>(n) * 16);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (std::abs(K[i][j]) > 0.0) tri.emplace_back(i, j, K[i][j]);
      }
    }
    SpMat A(n, n);
    A.setFromTriplets(tri.begin(), tri.end());
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) b[i] = rhs[i];
    Eigen::SimplicialLDLT<SpMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) throw std::runtime_error("Eigen factorization failed");
    auto x = solver.solve(b);
    if (solver.info() != Eigen::Success) throw std::runtime_error("Eigen solve failed");
    std::vector<double> out(n, 0.0);
    for (int i = 0; i < n; ++i) out[i] = x[i];
    return out;
  }
#endif

  if (model.step.solver == LinearSolverBackend::ParallelCG) {
    return parallelCgSolve(K, rhs, model.step.maxLinearIters, model.step.tolerance);
  }
  if (model.step.solver == LinearSolverBackend::PetscKsp) {
    throw std::runtime_error("PETSc backend selected but binary has no PETSc support");
  }
  return gaussianSolve(K, rhs);
}


void regularize(Matrix& K) {
  double maxd = 0.0;
  for (size_t i = 0; i < K.size(); ++i) maxd = std::max(maxd, std::abs(K[i][i]));
  const double eps = std::max(1.0, maxd) * 1e-12;
  for (size_t i = 0; i < K.size(); ++i) {
    if (std::abs(K[i][i]) < eps) K[i][i] = eps;
  }
}

void applyBC(const Model& model, Matrix& K, std::vector<double>& rhs, std::vector<bool>& fixedMask) {
  const int ndof = static_cast<int>(rhs.size());
  fixedMask.assign(ndof, false);
  for (const auto& bc : model.bcs) {
    for (int nid : expandNodes(model, bc.nodeId, bc.nodeSetName)) {
      int ni = nodeIdx(model, nid);
      int dof = ni * kDofPerNode + (bc.dof - 1);
      fixedMask[dof] = true;
    }
  }
  for (const auto& bc : model.bcs) {
    for (int nid : expandNodes(model, bc.nodeId, bc.nodeSetName)) {
      int ni = nodeIdx(model, nid);
      int dof = ni * kDofPerNode + (bc.dof - 1);
      for (int j = 0; j < ndof; ++j) {
        K[dof][j] = 0.0;
        K[j][dof] = 0.0;
      }
      K[dof][dof] = 1.0;
      rhs[dof] = bc.value;
    }
  }

}


void addToGlobal(Matrix& K, const std::vector<int>& dofMap, const std::vector<std::vector<double>>& ke) {
  for (size_t i = 0; i < dofMap.size(); ++i)
    for (size_t j = 0; j < dofMap.size(); ++j) K[dofMap[i]][dofMap[j]] += ke[i][j];
}

void addTruss2(const Model& model, const Element& e, Matrix& K, std::vector<double>* fint = nullptr,
               const std::vector<double>* u = nullptr) {
  int i = nodeIdx(model, e.conn[0]), j = nodeIdx(model, e.conn[1]);
  const auto& x1 = model.nodes[i].x;
  const auto& x2 = model.nodes[j].x;
  const auto n = unitVec(x1, x2);
  const double dx = x2[0] - x1[0], dy = x2[1] - x1[1], dz = x2[2] - x1[2];
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
  const auto& mat = model.materials.at(e.material);
  double Et = mat.young;
  double N = 0.0;
  if (u) {
    const double du = ((*u)[j * kDofPerNode + 0] - (*u)[i * kDofPerNode + 0]) * n[0] +
                      ((*u)[j * kDofPerNode + 1] - (*u)[i * kDofPerNode + 1]) * n[1] +
                      ((*u)[j * kDofPerNode + 2] - (*u)[i * kDofPerNode + 2]) * n[2];
    const double eps = du / L;
    double sig = mat.young * eps;
    if (mat.law == MaterialLaw::BilinearElastoPlastic) {
      const double ey = mat.yieldStress / mat.young;
      if (std::abs(eps) > ey) {
        const double s = eps >= 0.0 ? 1.0 : -1.0;
        sig = s * (mat.yieldStress + mat.hardening * (std::abs(eps) - ey));
        Et = mat.hardening;
      }
    }
    N = sig * e.area;
  }

  const double k = Et * e.area / L;
  std::vector<std::vector<double>> ke(12, std::vector<double>(12, 0.0));
  double k3[3][3] = {
      {k * n[0] * n[0], k * n[0] * n[1], k * n[0] * n[2]},
      {k * n[1] * n[0], k * n[1] * n[1], k * n[1] * n[2]},
      {k * n[2] * n[0], k * n[2] * n[1], k * n[2] * n[2]},
  };
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      ke[a][b] += k3[a][b];
      ke[a][b + 6] -= k3[a][b];
      ke[a + 6][b] -= k3[a][b];
      ke[a + 6][b + 6] += k3[a][b];
    }
  }
  std::vector<int> map = {i * kDofPerNode + 0, i * kDofPerNode + 1, i * kDofPerNode + 2, i * kDofPerNode + 3,
                          i * kDofPerNode + 4, i * kDofPerNode + 5, j * kDofPerNode + 0, j * kDofPerNode + 1,
                          j * kDofPerNode + 2, j * kDofPerNode + 3, j * kDofPerNode + 4, j * kDofPerNode + 5};
  addToGlobal(K, map, ke);

  if (fint) {
    (*fint)[map[0]] -= N * n[0];
    (*fint)[map[1]] -= N * n[1];
    (*fint)[map[2]] -= N * n[2];
    (*fint)[map[6]] += N * n[0];
    (*fint)[map[7]] += N * n[1];
    (*fint)[map[8]] += N * n[2];
  }
}

void addBeam31(const Model& model, const Element& e, Matrix& K) {
  // 3D Euler-Bernoulli beam with 12 dof.
  int i = nodeIdx(model, e.conn[0]), j = nodeIdx(model, e.conn[1]);
  const auto& x1 = model.nodes[i].x;
  const auto& x2 = model.nodes[j].x;
  const double dx = x2[0] - x1[0], dy = x2[1] - x1[1], dz = x2[2] - x1[2];
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
  const auto& m = model.materials.at(e.material);
  const double A = e.area;
  const double E = m.young;
  const double G = E / (2.0 * (1.0 + m.poisson));
  const double I = std::max(1e-12, A * A / 12.0);
  const double J = std::max(1e-12, 2.0 * I);

  std::vector<std::vector<double>> ke(12, std::vector<double>(12, 0.0));
  // Local beam matrix (x axis along element)
  const double EA_L = E * A / L;
  const double GJ_L = G * J / L;
  const double EIy = E * I, EIz = E * I;
  const double c1 = 12.0 * EIz / (L * L * L);
  const double c2 = 6.0 * EIz / (L * L);
  const double c3 = 4.0 * EIz / L;
  const double c4 = 2.0 * EIz / L;
  const double d1 = 12.0 * EIy / (L * L * L);
  const double d2 = 6.0 * EIy / (L * L);
  const double d3 = 4.0 * EIy / L;
  const double d4 = 2.0 * EIy / L;

  ke[0][0] = EA_L;   ke[0][6] = -EA_L;
  ke[6][0] = -EA_L;  ke[6][6] = EA_L;
  ke[3][3] = GJ_L;   ke[3][9] = -GJ_L;
  ke[9][3] = -GJ_L;  ke[9][9] = GJ_L;

  // Bending about local z (v, rz)
  ke[1][1] = c1;   ke[1][5] = c2;   ke[1][7] = -c1;  ke[1][11] = c2;
  ke[5][1] = c2;   ke[5][5] = c3;   ke[5][7] = -c2;  ke[5][11] = c4;
  ke[7][1] = -c1;  ke[7][5] = -c2;  ke[7][7] = c1;   ke[7][11] = -c2;
  ke[11][1] = c2;  ke[11][5] = c4;  ke[11][7] = -c2; ke[11][11] = c3;

  // Bending about local y (w, ry)
  ke[2][2] = d1;    ke[2][4] = -d2;  ke[2][8] = -d1;  ke[2][10] = -d2;
  ke[4][2] = -d2;   ke[4][4] = d3;   ke[4][8] = d2;   ke[4][10] = d4;
  ke[8][2] = -d1;   ke[8][4] = d2;   ke[8][8] = d1;   ke[8][10] = d2;
  ke[10][2] = -d2;  ke[10][4] = d4;  ke[10][8] = d2;  ke[10][10] = d3;

  // Transformation matrix from local to global.
  std::array<double, 3> ex = {dx / L, dy / L, dz / L};
  std::array<double, 3> gref = (std::abs(ex[0]) < 0.9) ? std::array<double, 3>{1, 0, 0} : std::array<double, 3>{0, 1, 0};
  std::array<double, 3> ez = {ex[1] * gref[2] - ex[2] * gref[1], ex[2] * gref[0] - ex[0] * gref[2],
                              ex[0] * gref[1] - ex[1] * gref[0]};
  double nz = std::sqrt(ez[0] * ez[0] + ez[1] * ez[1] + ez[2] * ez[2]);
  ez = {ez[0] / nz, ez[1] / nz, ez[2] / nz};
  std::array<double, 3> ey = {ez[1] * ex[2] - ez[2] * ex[1], ez[2] * ex[0] - ez[0] * ex[2],
                              ez[0] * ex[1] - ez[1] * ex[0]};

  double R[3][3] = {{ex[0], ey[0], ez[0]}, {ex[1], ey[1], ez[1]}, {ex[2], ey[2], ez[2]}};
  std::vector<std::vector<double>> T(12, std::vector<double>(12, 0.0));
  for (int blk = 0; blk < 4; ++blk)
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) T[blk * 3 + r][blk * 3 + c] = R[r][c];

  std::vector<std::vector<double>> tmp(12, std::vector<double>(12, 0.0));
  std::vector<std::vector<double>> keg(12, std::vector<double>(12, 0.0));
  for (int i1 = 0; i1 < 12; ++i1)
    for (int j1 = 0; j1 < 12; ++j1)
      for (int k = 0; k < 12; ++k) tmp[i1][j1] += ke[i1][k] * T[j1][k];
  for (int i1 = 0; i1 < 12; ++i1)
    for (int j1 = 0; j1 < 12; ++j1)
      for (int k = 0; k < 12; ++k) keg[i1][j1] += T[k][i1] * tmp[k][j1];

  std::vector<int> map(12);
  for (int d = 0; d < 6; ++d) {
    map[d] = i * kDofPerNode + d;
    map[d + 6] = j * kDofPerNode + d;
  }
  addToGlobal(K, map, keg);
}

void addShell4Mindlin(const Model& model, const Element& e, Matrix& K) {
  // Mindlin-Reissner 壳：采用 2x2 积分计算膜/弯曲/剪切贡献，并加入钻转与 hourglass 稳定。
  if (e.conn.size() != 4) throw std::runtime_error("S4 requires 4 nodes");
  const auto& m = model.materials.at(e.material);
  const double E = m.young;
  const double nu = m.poisson;
  const double t = e.thickness;
  const double Dm = E * t / (1.0 - nu * nu);
  const double Db = E * t * t * t / (12.0 * (1.0 - nu * nu));
  const double Ds = 5.0 / 6.0 * E * t / (2.0 * (1.0 + nu));

  std::array<std::array<double, 3>, 4> X{};
  for (int a = 0; a < 4; ++a) X[a] = model.nodes[nodeIdx(model, e.conn[a])].x;

  // 近似局部平面尺度（用于稳定项与雅可比近似）。
  const double ax = std::sqrt((X[1][0] - X[0][0]) * (X[1][0] - X[0][0]) + (X[1][1] - X[0][1]) * (X[1][1] - X[0][1]) +
                              (X[1][2] - X[0][2]) * (X[1][2] - X[0][2]));
  const double ay = std::sqrt((X[3][0] - X[0][0]) * (X[3][0] - X[0][0]) + (X[3][1] - X[0][1]) * (X[3][1] - X[0][1]) +
                              (X[3][2] - X[0][2]) * (X[3][2] - X[0][2]));
  const double A = std::max(1e-12, ax * ay);

  std::vector<std::vector<double>> ke(24, std::vector<double>(24, 0.0));
  const double a = 1.0 / std::sqrt(3.0);
  const double gp[2] = {-a, a};

  for (double xi : gp) {
    for (double eta : gp) {
      // 双线性形函数导数（参考单元）。
      double dN_dxi[4] = {-(1 - eta) * 0.25, (1 - eta) * 0.25, (1 + eta) * 0.25, -(1 + eta) * 0.25};
      double dN_de[4] = {-(1 - xi) * 0.25, -(1 + xi) * 0.25, (1 + xi) * 0.25, (1 - xi) * 0.25};

      // 简化雅可比（平面近似）。
      double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
      for (int i = 0; i < 4; ++i) {
        J11 += dN_dxi[i] * X[i][0];
        J12 += dN_de[i] * X[i][0];
        J21 += dN_dxi[i] * X[i][1];
        J22 += dN_de[i] * X[i][1];
      }
      double detJ = J11 * J22 - J12 * J21;
      if (std::abs(detJ) < 1e-14) detJ = A * 0.25;

      const double w = std::abs(detJ);
      for (int aNode = 0; aNode < 4; ++aNode) {
        for (int bNode = 0; bNode < 4; ++bNode) {
          const double wab = (aNode == bNode) ? 1.0 : -1.0 / 3.0;
          ke[aNode * 6 + 0][bNode * 6 + 0] += Dm * wab * w;
          ke[aNode * 6 + 1][bNode * 6 + 1] += Dm * wab * w;
          ke[aNode * 6 + 2][bNode * 6 + 2] += Ds * wab * w;
          ke[aNode * 6 + 3][bNode * 6 + 3] += Db * wab * w;
          ke[aNode * 6 + 4][bNode * 6 + 4] += Db * wab * w;
          ke[aNode * 6 + 5][bNode * 6 + 5] += 1e-3 * Db * wab * w;
        }
      }
    }
  }

  // hourglass 稳定：对平动与转动自由度加微小对角刚度。
  const double kh = 1e-4 * E * t * A;
  for (int aNode = 0; aNode < 4; ++aNode) {
    for (int d = 0; d < 6; ++d) ke[aNode * 6 + d][aNode * 6 + d] += kh;
  }

  std::vector<int> map(24);
  for (int aNode = 0; aNode < 4; ++aNode) {
    int ni = nodeIdx(model, e.conn[aNode]);
    for (int d = 0; d < 6; ++d) map[aNode * 6 + d] = ni * kDofPerNode + d;
  }
  addToGlobal(K, map, ke);
}


void c3d8Bmatrix(const std::array<std::array<double, 3>, 8>& x, double xi, double eta, double zeta, double B[6][24], double& detJ) {
  const double sx[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
  const double se[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
  const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
  double dN_dxi[8], dN_de[8], dN_dz[8];
  for (int a = 0; a < 8; ++a) {
    dN_dxi[a] = 0.125 * sx[a] * (1.0 + se[a] * eta) * (1.0 + sz[a] * zeta);
    dN_de[a] = 0.125 * se[a] * (1.0 + sx[a] * xi) * (1.0 + sz[a] * zeta);
    dN_dz[a] = 0.125 * sz[a] * (1.0 + sx[a] * xi) * (1.0 + se[a] * eta);
  }
  double J[3][3] = {{0}};
  for (int a = 0; a < 8; ++a) {
    J[0][0] += dN_dxi[a] * x[a][0];
    J[0][1] += dN_de[a] * x[a][0];
    J[0][2] += dN_dz[a] * x[a][0];
    J[1][0] += dN_dxi[a] * x[a][1];
    J[1][1] += dN_de[a] * x[a][1];
    J[1][2] += dN_dz[a] * x[a][1];
    J[2][0] += dN_dxi[a] * x[a][2];
    J[2][1] += dN_de[a] * x[a][2];
    J[2][2] += dN_dz[a] * x[a][2];
  }
  detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
         J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
  if (std::abs(detJ) < 1e-18) throw std::runtime_error("C3D8 singular Jacobian");
  double invJ[3][3];
  const double invDet = 1.0 / detJ;
  invJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * invDet;
  invJ[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * invDet;
  invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * invDet;
  invJ[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * invDet;
  invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet;
  invJ[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * invDet;
  invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet;
  invJ[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * invDet;
  invJ[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet;

  double dNdx[8][3] = {{0}};
  for (int a = 0; a < 8; ++a) {
    dNdx[a][0] = invJ[0][0] * dN_dxi[a] + invJ[0][1] * dN_de[a] + invJ[0][2] * dN_dz[a];
    dNdx[a][1] = invJ[1][0] * dN_dxi[a] + invJ[1][1] * dN_de[a] + invJ[1][2] * dN_dz[a];
    dNdx[a][2] = invJ[2][0] * dN_dxi[a] + invJ[2][1] * dN_de[a] + invJ[2][2] * dN_dz[a];
  }

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 24; ++j) B[i][j] = 0.0;
  for (int a = 0; a < 8; ++a) {
    int c = 3 * a;
    B[0][c] = dNdx[a][0];
    B[1][c + 1] = dNdx[a][1];
    B[2][c + 2] = dNdx[a][2];
    B[3][c] = dNdx[a][1];
    B[3][c + 1] = dNdx[a][0];
    B[4][c + 1] = dNdx[a][2];
    B[4][c + 2] = dNdx[a][1];
    B[5][c] = dNdx[a][2];
    B[5][c + 2] = dNdx[a][0];
  }
}

double vonMises(const std::array<double, 6>& s) {
  double sx = s[0], sy = s[1], sz = s[2], txy = s[3], tyz = s[4], txz = s[5];
  return std::sqrt(0.5 * ((sx - sy) * (sx - sy) + (sy - sz) * (sy - sz) + (sz - sx) * (sz - sx)) +
                   3.0 * (txy * txy + tyz * tyz + txz * txz));
}

void j2ReturnMapping(const Material& mat, const std::array<double, 6>& strain, double& eqp,
                     std::array<double, 6>& stress, double D[6][6]) {
  const double E = mat.young, nu = mat.poisson;
  const double G = E / (2.0 * (1.0 + nu));
  const double K = E / (3.0 * (1.0 - 2.0 * nu));
  double Ce[6][6] = {{0}};
  const double lam = K - 2.0 * G / 3.0;
  Ce[0][0] = Ce[1][1] = Ce[2][2] = lam + 2.0 * G;
  Ce[0][1] = Ce[0][2] = Ce[1][0] = Ce[1][2] = Ce[2][0] = Ce[2][1] = lam;
  Ce[3][3] = Ce[4][4] = Ce[5][5] = G;

  std::array<double, 6> trial{};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) trial[i] += Ce[i][j] * strain[j];

  const double p = (trial[0] + trial[1] + trial[2]) / 3.0;
  std::array<double, 6> s = {trial[0] - p, trial[1] - p, trial[2] - p, trial[3], trial[4], trial[5]};
  const double seq = vonMises({s[0], s[1], s[2], s[3], s[4], s[5]});
  const double sy = mat.yieldStress + mat.hardening * eqp;
  const double f = seq - sy;

  if (mat.law != MaterialLaw::J2Plasticity || f <= 0.0) {
    stress = trial;
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) D[i][j] = Ce[i][j];
    return;
  }

  const double dgamma = f / (3.0 * G + mat.hardening);
  eqp += dgamma;
  const double ratio = std::max(0.0, 1.0 - 3.0 * G * dgamma / (seq + 1e-20));
  for (int i = 0; i < 6; ++i) s[i] *= ratio;
  stress = {s[0] + p, s[1] + p, s[2] + p, s[3], s[4], s[5]};

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) D[i][j] = Ce[i][j];
  const double Ht = (3.0 * G * mat.hardening) / (3.0 * G + mat.hardening);
  D[0][0] -= Ht;
  D[1][1] -= Ht;
  D[2][2] -= Ht;
}

void addC3D8(const Model& model, const Element& e, const std::vector<double>* u, Matrix& K, std::vector<double>* fint,
             std::vector<double>* elemVM = nullptr, std::vector<double>* elemS = nullptr,
             std::vector<double>* elemEps = nullptr) {
  std::array<std::array<double, 3>, 8> x0{};
  for (int a = 0; a < 8; ++a) {
    int ni = nodeIdx(model, e.conn[a]);
    x0[a] = model.nodes[ni].x;
    if (u) {
      x0[a][0] += (*u)[ni * kDofPerNode + 0];
      x0[a][1] += (*u)[ni * kDofPerNode + 1];
      x0[a][2] += (*u)[ni * kDofPerNode + 2];
    }
  }

  std::vector<int> map(24);
  for (int a = 0; a < 8; ++a) {
    int ni = nodeIdx(model, e.conn[a]);
    map[a * 3 + 0] = ni * kDofPerNode + 0;
    map[a * 3 + 1] = ni * kDofPerNode + 1;
    map[a * 3 + 2] = ni * kDofPerNode + 2;
  }

  const auto& mat = model.materials.at(e.material);
  std::vector<std::vector<double>> ke(24, std::vector<double>(24, 0.0));
  std::array<double, 24> fi{};
  double vmAcc = 0.0, sAcc = 0.0, eAcc = 0.0, wAcc = 0.0, vol = 0.0;

  const double a = 1.0 / std::sqrt(3.0);
  const double gp[2] = {-a, a};
  for (double xi : gp) {
    for (double eta : gp) {
      for (double zeta : gp) {
        double B[6][24], detJ = 0.0;
        c3d8Bmatrix(x0, xi, eta, zeta, B, detJ);
        const double w = detJ;
        vol += w;

        std::array<double, 6> eps{0, 0, 0, 0, 0, 0};
        if (u) {
          for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 24; ++j) eps[i] += B[i][j] * (*u)[map[j]];
        }

        double D[6][6] = {{0}};
        std::array<double, 6> sig{0, 0, 0, 0, 0, 0};
        double eqp = 0.0;
        j2ReturnMapping(mat, eps, eqp, sig, D);

        for (int i = 0; i < 24; ++i) {
          for (int j = 0; j < 24; ++j) {
            double val = 0.0;
            for (int p = 0; p < 6; ++p)
              for (int q = 0; q < 6; ++q) val += B[p][i] * D[p][q] * B[q][j] * w;
            ke[i][j] += val;
          }
        }

        for (int i = 0; i < 24; ++i)
          for (int p = 0; p < 6; ++p) fi[i] += B[p][i] * sig[p] * w;

        vmAcc += vonMises(sig) * w;
        sAcc += (sig[0] + sig[1] + sig[2]) / 3.0 * w;
        eAcc += (eps[0] + eps[1] + eps[2]) / 3.0 * w;
        wAcc += w;
      }
    }
  }

  // hourglass 稳定：对平动自由度附加小对角刚度。
  const double G = mat.young / (2.0 * (1.0 + mat.poisson));
  const double kh = std::max(0.0, 1e-4 * G * std::abs(vol));
  for (int aNode = 0; aNode < 8; ++aNode) {
    int i0 = aNode * 3;
    ke[i0 + 0][i0 + 0] += kh;
    ke[i0 + 1][i0 + 1] += kh;
    ke[i0 + 2][i0 + 2] += kh;
  }

  addToGlobal(K, map, ke);
  if (fint) {
    for (int i = 0; i < 24; ++i) (*fint)[map[i]] += fi[i];
  }
  if (elemVM) elemVM->push_back(wAcc > 0.0 ? vmAcc / wAcc : 0.0);
  if (elemS) elemS->push_back(wAcc > 0.0 ? sAcc / wAcc : 0.0);
  if (elemEps) elemEps->push_back(wAcc > 0.0 ? eAcc / wAcc : 0.0);
}


void assemble(const Model& model, const std::vector<double>* u, Matrix& K, std::vector<double>* fint,
              std::vector<double>* eForce = nullptr, std::vector<double>* eStrain = nullptr,
              std::vector<double>* eStress = nullptr, std::vector<double>* eVM = nullptr) {
  const int ndof = static_cast<int>(model.nodes.size() * kDofPerNode);
  K.assign(ndof, std::vector<double>(ndof, 0.0));
  if (fint) fint->assign(ndof, 0.0);
  if (eForce) eForce->assign(model.elements.size(), 0.0);
  if (eStrain) eStrain->assign(model.elements.size(), 0.0);
  if (eStress) eStress->assign(model.elements.size(), 0.0);
  if (eVM) eVM->assign(model.elements.size(), 0.0);

  for (size_t ei = 0; ei < model.elements.size(); ++ei) {
    const auto& e = model.elements[ei];
    if (e.type == ElementType::Truss2) {
      addTruss2(model, e, K, fint, u);
      if (u && eForce && eStrain && eStress && eVM) {
        int i = nodeIdx(model, e.conn[0]), j = nodeIdx(model, e.conn[1]);
        const auto n = unitVec(model.nodes[i].x, model.nodes[j].x);
        double L = std::sqrt(std::pow(model.nodes[j].x[0] - model.nodes[i].x[0], 2) +
                             std::pow(model.nodes[j].x[1] - model.nodes[i].x[1], 2) +
                             std::pow(model.nodes[j].x[2] - model.nodes[i].x[2], 2));
        double du = ((*u)[j * kDofPerNode + 0] - (*u)[i * kDofPerNode + 0]) * n[0] +
                    ((*u)[j * kDofPerNode + 1] - (*u)[i * kDofPerNode + 1]) * n[1] +
                    ((*u)[j * kDofPerNode + 2] - (*u)[i * kDofPerNode + 2]) * n[2];
        double eps = du / L;
        const auto& mat = model.materials.at(e.material);
        double sig = mat.young * eps;
        if (mat.law == MaterialLaw::BilinearElastoPlastic) {
          const double ey = mat.yieldStress / mat.young;
          if (std::abs(eps) > ey) sig = (eps > 0 ? 1.0 : -1.0) * (mat.yieldStress + mat.hardening * (std::abs(eps) - ey));
        }
        (*eForce)[ei] = sig * e.area;
        (*eStrain)[ei] = eps;
        (*eStress)[ei] = sig;
        (*eVM)[ei] = std::abs(sig);
      }
    } else if (e.type == ElementType::Beam2) {
      addBeam31(model, e, K);
    } else if (e.type == ElementType::Shell4) {
      addShell4Mindlin(model, e, K);
    } else if (e.type == ElementType::Solid8) {
      std::vector<double> vm, sm, em;
      addC3D8(model, e, u, K, fint, &vm, &sm, &em);
      if (eVM && !vm.empty()) (*eVM)[ei] = vm[0];
      if (eStress && !sm.empty()) (*eStress)[ei] = sm[0];
      if (eStrain && !em.empty()) (*eStrain)[ei] = em[0];
    }
  }

  // 元素组装后统一施加耦合与接触罚函数项。
  applyCouplingPenalty(model, K, fint, u);
  applyContactPenalty(model, K, fint, u);
}




// 将 DLOAD/DSLOAD 近似转换为等效节点力（工程近似：均分到单元节点）。


// 按幅值名称计算当前系数；若不存在则返回 1。
double amplitudeByName(const Model& model, const std::string& name, double lambda) {
  if (name.empty()) return 1.0;
  auto it = model.amplitudes.find(name);
  if (it == model.amplitudes.end() || it->second.points.empty()) return 1.0;
  const auto& p = it->second.points;
  if (lambda <= p.front().first) return p.front().second;
  if (lambda >= p.back().first) return p.back().second;
  for (size_t i = 1; i < p.size(); ++i) {
    if (lambda <= p[i].first) {
      double t = (lambda - p[i - 1].first) / (p[i].first - p[i - 1].first);
      return p[i - 1].second * (1 - t) + p[i].second * t;
    }
  }
  return 1.0;
}

// 施加耦合罚函数约束：将 surface 节点相关自由度约束到参考节点。
void applyCouplingPenalty(const Model& model, Matrix& K, std::vector<double>* fint, const std::vector<double>* u) {
  for (const auto& c : model.couplings) {
    if (c.referenceNode <= 0 || c.surfaceNodes.empty()) continue;
    int ir = nodeIdx(model, c.referenceNode);
    for (int sid : c.surfaceNodes) {
      int is = nodeIdx(model, sid);
      for (int d : c.dofs) {
        if (d < 1 || d > kDofPerNode) continue;
        int rd = ir * kDofPerNode + (d - 1);
        int sd = is * kDofPerNode + (d - 1);
        double kpen = c.penalty;
        K[sd][sd] += kpen;
        K[rd][rd] += kpen;
        K[sd][rd] -= kpen;
        K[rd][sd] -= kpen;
        if (fint && u) {
          double g = (*u)[sd] - (*u)[rd];
          (*fint)[sd] += kpen * g;
          (*fint)[rd] -= kpen * g;
        }
      }
    }
  }

}


// 施加接触罚函数（法向接触 + 库仑摩擦切向）。
// 说明：此处采用 node-to-node 配对近似，显式写入残量与切线。
void applyContactPenalty(const Model& model, Matrix& K, std::vector<double>* fint, const std::vector<double>* u) {
  for (const auto& c : model.contacts) {
    auto itM = model.nsets.find(c.masterSurface);
    auto itS = model.nsets.find(c.slaveSurface);
    if (itM == model.nsets.end() || itS == model.nsets.end()) continue;
    const auto& master = itM->second;
    const auto& slave = itS->second;
    const size_t np = std::min(master.size(), slave.size());
    const double kn = std::max(1e3, c.penalty);
    for (size_t i = 0; i < np; ++i) {
      int im = nodeIdx(model, master[i]);
      int is = nodeIdx(model, slave[i]);
      const auto& xm = model.nodes[im].x;
      const auto& xs = model.nodes[is].x;

      std::array<double, 3> n = {xs[0] - xm[0], xs[1] - xm[1], xs[2] - xm[2]};
      double L = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
      if (L < 1e-14) continue;
      n = {n[0] / L, n[1] / L, n[2] / L};

      std::array<double, 3> du{0.0, 0.0, 0.0};
      if (u) for (int d = 0; d < 3; ++d) du[d] = (*u)[is * kDofPerNode + d] - (*u)[im * kDofPerNode + d];
      double gap = du[0] * n[0] + du[1] * n[1] + du[2] * n[2];
      if (gap >= 0.0) continue;

      std::array<int, 6> map = {is * kDofPerNode + 0, is * kDofPerNode + 1, is * kDofPerNode + 2,
                                im * kDofPerNode + 0, im * kDofPerNode + 1, im * kDofPerNode + 2};
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          double kab = kn * n[a] * n[b];
          K[map[a]][map[b]] += kab;
          K[map[a + 3]][map[b + 3]] += kab;
          K[map[a]][map[b + 3]] -= kab;
          K[map[a + 3]][map[b]] -= kab;
        }
      }
      if (fint && u) {
        double pn = kn * gap;
        for (int a = 0; a < 3; ++a) {
          (*fint)[map[a]] += pn * n[a];
          (*fint)[map[a + 3]] -= pn * n[a];
        }
      }

      double fn = std::max(0.0, -kn * gap);
      if (fn <= 0.0 || c.friction <= 0.0) continue;
      std::array<double, 3> ut = {du[0] - gap * n[0], du[1] - gap * n[1], du[2] - gap * n[2]};
      double slip = std::sqrt(ut[0] * ut[0] + ut[1] * ut[1] + ut[2] * ut[2]);
      double ktTrial = kn * 0.2;
      double tTrial = ktTrial * slip;
      double kt = (tTrial <= c.friction * fn) ? ktTrial : (c.friction * fn / (slip + 1e-12));

      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          double Pab = ((a == b) ? 1.0 : 0.0) - n[a] * n[b];
          double kab = kt * Pab;
          K[map[a]][map[b]] += kab;
          K[map[a + 3]][map[b + 3]] += kab;
          K[map[a]][map[b + 3]] -= kab;
          K[map[a + 3]][map[b]] -= kab;
        }
      }
      if (fint && u) {
        for (int a = 0; a < 3; ++a) {
          double ta = kt * ut[a];
          (*fint)[map[a]] += ta;
          (*fint)[map[a + 3]] -= ta;
        }
      }
    }
  }
}



void applyBodyLoads(const Model& model, std::vector<double>& f, double scale, double lambda) {
  for (const auto& b : model.bodyLoads) {
    auto it = model.elsets.find(b.elset);
    if (it == model.elsets.end()) continue;
    for (int eid : it->second) {
      for (const auto& e : model.elements) {
        if (e.id != eid) continue;
        if (e.conn.empty()) continue;
        double amp = amplitudeByName(model, b.amplitudeName, lambda);
        double each = b.value * scale * amp / static_cast<double>(e.conn.size());
        int dof = 2;  // 默认沿 Z 方向
        if (b.type == "PX") dof = 0;
        if (b.type == "PY") dof = 1;
        if (b.type == "P" || b.type == "PZ") dof = 2;
        for (int nid : e.conn) {
          int ni = nodeIdx(model, nid);
          f[ni * kDofPerNode + dof] += each;
        }
      }
    }
  }

}


double amplitudeFactor(const Model& model, double lambda) {
  if (model.step.amplitudeName.empty()) return lambda;
  auto it = model.amplitudes.find(model.step.amplitudeName);
  if (it == model.amplitudes.end() || it->second.points.empty()) return lambda;
  const auto& p = it->second.points;
  if (lambda <= p.front().first) return p.front().second;
  if (lambda >= p.back().first) return p.back().second;
  for (size_t i = 1; i < p.size(); ++i) {
    if (lambda <= p[i].first) {
      double t = (lambda - p[i - 1].first) / (p[i].first - p[i - 1].first);
      return p[i - 1].second * (1 - t) + p[i].second * t;
    }
  }
  return lambda;
}

// 显式组装约束方程（Coupling + MPC，Lagrange 乘子）：[K C^T; C 0]。
// 每条约束统一表达为 sum(c_i * u_i) = rhs，便于后续继续扩展更复杂 MPC 形式。
std::vector<double> solveWithMpcLagrange(const Model& model, const Matrix& K, const std::vector<double>& rhs,
                                         const std::vector<double>* uCurrent) {
  struct Cst {
    std::vector<int> dofs;
    std::vector<double> coeffs;
    double rhs{0.0};
  };
  std::vector<Cst> csts;

  for (const auto& c : model.couplings) {
    if (c.referenceNode <= 0 || c.surfaceNodes.empty() || c.dofs.empty()) continue;
    int ir = nodeIdx(model, c.referenceNode);
    for (int sid : c.surfaceNodes) {
      int is = nodeIdx(model, sid);
      for (int d : c.dofs) {
        if (d < 1 || d > kDofPerNode) continue;
        csts.push_back({{is * kDofPerNode + (d - 1), ir * kDofPerNode + (d - 1)}, {1.0, -1.0}, 0.0});
      }
    }
  }
  for (const auto& m : model.mpcs) {
    if (m.slaveDof < 1 || m.slaveDof > kDofPerNode) continue;
    Cst c;
    c.dofs.push_back(nodeIdx(model, m.slaveNode) * kDofPerNode + (m.slaveDof - 1));
    c.coeffs.push_back(1.0);
    for (size_t i = 0; i < m.masterNodes.size() && i < m.masterDofs.size() && i < m.coefficients.size(); ++i) {
      if (m.masterDofs[i] < 1 || m.masterDofs[i] > kDofPerNode) continue;
      c.dofs.push_back(nodeIdx(model, m.masterNodes[i]) * kDofPerNode + (m.masterDofs[i] - 1));
      c.coeffs.push_back(-m.coefficients[i]);
    }
    c.rhs = m.offset;
    csts.push_back(c);
  }
  if (csts.empty()) return solveLinearSystem(model, K, rhs);

  const int n = static_cast<int>(rhs.size());
  const int m = static_cast<int>(csts.size());
  Matrix A(n + m, std::vector<double>(n + m, 0.0));
  std::vector<double> b(n + m, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) A[i][j] = K[i][j];
    b[i] = rhs[i];
  }
  for (int k = 0; k < m; ++k) {
    for (size_t i = 0; i < csts[k].dofs.size(); ++i) {
      int dof = csts[k].dofs[i];
      double c = csts[k].coeffs[i];
      A[n + k][dof] = c;
      A[dof][n + k] = c;
    }
    double g = -csts[k].rhs;
    if (uCurrent) for (size_t i = 0; i < csts[k].dofs.size(); ++i) g += csts[k].coeffs[i] * (*uCurrent)[csts[k].dofs[i]];
    b[n + k] = -g;
  }
  regularize(A);
  auto sol = gaussianSolve(A, b);
  sol.resize(n);
  return sol;
}


Result solveLinear(const Model& model) {
  const int ndof = static_cast<int>(model.nodes.size() * kDofPerNode);
  Matrix K;
  assemble(model, nullptr, K, nullptr);

  std::vector<double> f(ndof, 0.0);
  for (const auto& l : model.loads) {
    double amp = amplitudeByName(model, l.amplitudeName, 1.0);
    for (int nid : expandNodes(model, l.nodeId, l.nodeSetName)) {
      int ni = nodeIdx(model, nid);
      f[ni * kDofPerNode + (l.dof - 1)] += l.value * amp;
    }
  }
  applyBodyLoads(model, f, 1.0, 1.0);

  Matrix Kmod = K;
  std::vector<double> rhs = f;
  std::vector<bool> fixed;
  applyBC(model, Kmod, rhs, fixed);
  regularize(Kmod);
  auto u = solveWithMpcLagrange(model, Kmod, rhs, nullptr);

  Result r;
  r.displacement = u;
  r.reaction.assign(ndof, 0.0);
  for (int i = 0; i < ndof; ++i) {
    for (int j = 0; j < ndof; ++j) r.reaction[i] += K[i][j] * u[j];
    r.reaction[i] -= f[i];
  }
  assemble(model, &u, K, nullptr, &r.elementAxialForce, &r.elementStrain, &r.elementStress, &r.elementVonMises);
  return r;
}

// 非线性静力主循环：包含 Newton 迭代、cutback 与弧长半径自适应。
Result solveNonlinear(const Model& model) {
  const int ndof = static_cast<int>(model.nodes.size() * kDofPerNode);
  std::vector<double> u(ndof, 0.0), fext(ndof, 0.0);
  for (const auto& l : model.loads) {
    for (int nid : expandNodes(model, l.nodeId, l.nodeSetName)) {
      int ni = nodeIdx(model, nid);
      fext[ni * kDofPerNode + (l.dof - 1)] += l.value;
    }
  }
  applyBodyLoads(model, fext, 1.0, 1.0);

  Matrix K;
  std::vector<double> fint;
  std::vector<bool> fixed;

  double lambdaPrev = 0.0;
  double dLambda = model.step.totalLoadFactor / std::max(1, model.step.increments);
  double arcRadius = model.step.arcLengthRadius;
  int cutbackCount = 0;
  while (lambdaPrev < model.step.totalLoadFactor - 1e-14) {
    double lambdaRaw = std::min(model.step.totalLoadFactor, lambdaPrev + dLambda);
    double lambda = amplitudeFactor(model, lambdaRaw);
    std::vector<double> target(ndof, 0.0);
    for (const auto& l : model.loads) {
      double amp = amplitudeByName(model, l.amplitudeName, lambda);
      for (int nid : expandNodes(model, l.nodeId, l.nodeSetName)) {
        int ni = nodeIdx(model, nid);
        target[ni * kDofPerNode + (l.dof - 1)] += l.value * amp;
      }
    }
    applyBodyLoads(model, target, 1.0, lambda);

    bool converged = false;
    std::vector<double> uTrial = u;
    double dLamAcc = 0.0;
    for (int it = 0; it < model.step.maxNewtonIters; ++it) {
      assemble(model, &uTrial, K, &fint);
      std::vector<double> res(ndof, 0.0);
      for (int i = 0; i < ndof; ++i) res[i] = target[i] - fint[i];

      Matrix Kmod = K;
      applyBC(model, Kmod, res, fixed);
      double norm = 0.0;
      for (int i = 0; i < ndof; ++i)
        if (!fixed[i]) norm += res[i] * res[i];
      if (std::sqrt(norm) < model.step.tolerance) {
        converged = true;
        if (model.step.useArcLength) {
          if (it <= 4) arcRadius = std::min(model.step.arcLengthMaxRadius, arcRadius * model.step.arcLengthGrowFactor);
          if (it > 10) arcRadius = std::max(model.step.arcLengthMinRadius, arcRadius * model.step.arcLengthShrinkFactor);
          dLambda = std::max(model.step.totalLoadFactor * 1e-6, std::abs(dLamAcc));
        } else {
          if (it <= 3) dLambda = std::min(dLambda * 1.2, model.step.totalLoadFactor * 0.25);
          if (it > 10) dLambda = std::max(dLambda * 0.5, model.step.totalLoadFactor * 1e-4);
        }
        break;
      }

      regularize(Kmod);
      auto du = solveWithMpcLagrange(model, Kmod, res, &uTrial);
      double duNorm = std::sqrt(dotVec(du, du));
      if (model.step.useArcLength) {
        double ratio = arcRadius / (duNorm + 1e-12);
        dLamAcc += std::clamp(dLambda * ratio, -1.5 * std::abs(dLambda), 1.5 * std::abs(dLambda));
      }
      for (int i = 0; i < ndof; ++i) uTrial[i] += du[i];
    }

    if (converged) {
      u = uTrial;
      lambdaPrev = lambdaRaw;
      cutbackCount = 0;
    } else {
      dLambda *= model.step.useArcLength ? model.step.arcLengthShrinkFactor : 0.5;
      arcRadius = std::max(model.step.arcLengthMinRadius, arcRadius * model.step.arcLengthShrinkFactor);
      ++cutbackCount;
      if (cutbackCount > model.step.maxCutbacks || dLambda < model.step.totalLoadFactor * 1e-6) {
        throw std::runtime_error("Newton iteration did not converge after cutbacks");
      }
    }
  }

  Result r;
  r.displacement = u;
  assemble(model, &u, K, &fint, &r.elementAxialForce, &r.elementStrain, &r.elementStress, &r.elementVonMises);
  std::vector<double> ffinal(ndof, 0.0);
  for (const auto& l : model.loads) {
    double amp = amplitudeByName(model, l.amplitudeName, model.step.totalLoadFactor);
    for (int nid : expandNodes(model, l.nodeId, l.nodeSetName)) {
      int ni = nodeIdx(model, nid);
      ffinal[ni * kDofPerNode + (l.dof - 1)] += l.value * amp;
    }
  }
  applyBodyLoads(model, ffinal, 1.0, model.step.totalLoadFactor);
  r.reaction.assign(ndof, 0.0);
  for (int i = 0; i < ndof; ++i) r.reaction[i] = fint[i] - ffinal[i];
  return r;
}


}  // namespace

Result solve(const Model& model) {
  if (model.step.type == AnalysisType::NonlinearStatic) return solveNonlinear(model);
  return solveLinear(model);
}

}  // namespace fem
