# C++ 静力学有限元求解器（Abaqus inp -> VTK）

本轮在上一版基础上继续增强，目标进一步逼近 Abaqus 静力学流程，并保持工程可持续迭代。

## 1) 已增强内容（本次）

- **并行机制引入**
  - 新增 OpenMP 并行路径（矩阵向量乘、向量点积、并行 CG 迭代）。
  - 新增求解后端：`pcg`（Parallel Conjugate Gradient）。

- **载荷与边界表达增强**
  - `*BOUNDARY` / `*CLOAD` 同时支持“节点号”与“节点集名称”。
  - 新增 `*DSLOAD` 解析，并与 `*DLOAD` 统一进入体载荷容器。
  - 在求解器中增加 `DLOAD/DSLOAD` 到等效节点力的近似转换（工程均分策略）。

- **Abaqus 关键字继续补全（解析层）**
  - 继续支持并完善：`*Amplitude`、`*Coupling`、`*Kinematic`、`*Contact Pair`、`*Node Output`、`*Element Output`。

- **注释规范**
  - 本轮新增与改动代码均使用中文注释，便于团队维护。

## 2) 当前能力总览

- 自由度：每节点 6 自由度（UX, UY, UZ, RX, RY, RZ）
- 单元：
  - `T3D2/T2D2`：2 节点桁架
  - `B31/B33`：3D Euler-Bernoulli 梁（12 DOF）
  - `S4/S4R`：Mindlin-Reissner 壳（简化 4 节点 24 DOF）
  - `C3D8/C3D8R`：8 节点实体（中心点积分）
- 材料：
  - 线弹性
  - 双线性弹塑性（杆类）
  - J2 塑性（C3D8 单积分点 return-mapping）
- 分析：线性与非线性静力（Newton-Raphson）
- 线性求解后端：Dense / Eigen / PETSc(KSP) / ParallelCG(OpenMP)
- 输出：VTK（位移、反力、轴力、应变、应力、von_mises）

## 3) 关键字支持

### 已支持

- `*Node`
- `*Element, type=...`
- `*Nset`
- `*Elset`
- `*Material, name=...`
- `*Elastic`
- `*Plastic`
- `*Truss Section`
- `*Beam Section`
- `*Shell Section`
- `*Boundary`（节点号 / 节点集）
- `*Cload`（节点号 / 节点集）
- `*Dload`
- `*Dsload`
- `*Amplitude`
- `*Step`（含 `nlgeom=YES` 与 `amplitude=...`）
- `*Static`
- `*Controls`
- `*Coupling`
- `*Kinematic`
- `*Contact Pair`
- `*Node Output`
- `*Element Output`
- `*End Step`

### 当前说明

- 接触、耦合、输出请求已完成解析与模型存储。
- 接触/耦合完整残量与切线刚度（罚函数/拉格朗日乘子）可在现有模型层继续增强。

## 4) macOS 安装与使用（详细）

### 4.1 依赖安装

```bash
brew install cmake eigen
# 可选：PETSc
brew install petsc
```

### 4.2 编译

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=ON -DFEM_ENABLE_PETSC=ON -DFEM_ENABLE_OPENMP=ON
cmake --build build -j
```

### 4.3 运行示例

```bash
./build/fem_solver examples/linear_truss.inp output/linear.vtk 1.0 dense
./build/fem_solver examples/nonlinear_truss.inp output/nonlinear.vtk 1.0 pcg
./build/fem_solver examples/c3d8_linear.inp output/c3d8_linear.vtk 1.0 eigen
./build/fem_solver examples/truss_plastic.inp output/truss_plastic.vtk 1.0 petsc
```

## 5) Windows 安装与使用（详细）

### 5.1 推荐环境

- Visual Studio 2022（C++ 工作负载）
- CMake
- 可选：Eigen（vcpkg）、PETSc（MSYS2/自编译）

### 5.2 编译

```bat
cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=ON -DFEM_ENABLE_PETSC=OFF -DFEM_ENABLE_OPENMP=ON
cmake --build build -j
```

### 5.3 运行

```bat
build\fem_solver.exe examples\linear_truss.inp output\linear.vtk 1.0 dense
build\fem_solver.exe examples\advanced_keywords.inp output\advanced.vtk 1.0 pcg
```

## 6) 命令行

```bash
fem_solver <input.inp> <output.vtk> [deform_scale] [solver=dense|eigen|petsc|pcg]
```

## 7) 说明与后续建议

- 该版本已具备“并行 + 多单元 + 非线性 + 更完整 inp 解析”的可运行骨架。
- 若要进一步对标 Abaqus 静力学，建议下一步优先实现：
  1. 接触残量/切线刚度（法向接触 + 摩擦）
  2. 耦合约束方程显式组装（MPC/Lagrange）
  3. 壳与实体的多积分点与 hourglass 控制
  4. 非线性步自适应弧长法（支持 cutback + RIKS 风格半径调节）

本轮新增：
- `*Contact Pair` 支持 `friction=` 参数输入摩擦系数。
- 非线性求解新增 cutback 机制（增量失败自动减步重算）。
