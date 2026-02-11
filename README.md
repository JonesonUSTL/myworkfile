# mini_fem_solver（Abaqus `.inp` -> VTK）

一个面向教学与工程原型迭代的 C++ 静力学有限元求解器。目标是逐步对标 Abaqus 的非线性结构静力流程，并保持代码可读、可扩展、可跨平台部署。

---

## 1. 你能直接做什么

- 读取 Abaqus 风格 `.inp` 输入文件并求解。
- 输出 VTK 文件，用 ParaView 直接查看变形与场变量。
- 支持线性与非线性静力计算（Newton-Raphson + cutback + 弧长参数）。
- 支持接触（法向 + 摩擦罚函数）、Coupling/Kinematic 和显式 MPC 约束。
- 支持 macOS（含 Apple Silicon/M4）与 Windows 的一键部署脚本。

---

## 2. 快速开始（5 分钟）

### 2.1 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=ON -DFEM_ENABLE_PETSC=OFF -DFEM_ENABLE_OPENMP=ON
cmake --build build -j
```

### 2.2 查看帮助/示例

```bash
./build/fem_solver --help
./build/fem_solver --examples
./build/fem_solver --solvers
```

### 2.3 跑一个最小算例

```bash
./build/fem_solver examples/classic_bar_tension.inp output/classic_bar.vtk 1.0 dense
```

---

## 3. 命令行与交互

```bash
fem_solver <input.inp> <output.vtk> [deform_scale] [solver=dense|eigen|petsc|pcg]
```

### 3.1 功能命令

- `--help`：打印完整帮助。
- `--examples`：打印内置案例列表。
- `--version`：打印程序版本。
- `--solvers`：打印支持的线性求解器后端。
- `--check-inp <input.inp>`：仅解析输入并给出模型统计（不求解）。

### 3.2 典型用法

```bash
./build/fem_solver --check-inp examples/advanced_keywords.inp
./build/fem_solver examples/linear_truss.inp output/linear.vtk 1.0 dense
./build/fem_solver examples/nonlinear_truss.inp output/nonlinear.vtk 1.0 pcg
```

---

## 4. 已支持能力

### 4.1 单元与材料

- 单元：`T3D2/T2D2`、`B31/B33`、`S4/S4R`、`C3D8/C3D8R`
- 材料：线弹性、双线性弹塑性（杆类）、J2（实体）

### 4.2 关键字

- 网格：`*Node` `*Element` `*Nset` `*Elset`
- 材料/截面：`*Material` `*Elastic` `*Plastic` `*Truss Section` `*Beam Section` `*Shell Section`
- 载荷/边界：`*Boundary` `*Cload` `*Dload` `*Dsload` `*Amplitude`
- 分析步：`*Step` `*Static` `*Controls`
- 非线性增强：`*Static, riks`
- 约束/接触：`*Coupling` `*Kinematic` `*MPC` `*Contact Pair`
- 输出请求：`*Output, field, frequency=` `*Field Output` `*Node Output` `*Element Output` `*End Step`
- 模型组织：`*Heading`（模型名）`*Part` `*End Part`
- 材料扩展：`*Density` `*Expansion` `*Conductivity` `*Specific Heat`
- 分析步扩展：`*Static`（初始/总量/最小/最大增量）

### 4.3 求解策略

- 线性：Dense / Eigen / PETSc / OpenMP PCG
- 非线性：Newton-Raphson、cutback、弧长半径自适应参数
- 约束：Lagrange 增广方程（Coupling + MPC）
- 接触：法向罚函数 + 库仑摩擦切向投影

---

## 5. 经典案例与回归建议

### 5.1 入门三件套（推荐每次改动后先跑）

```bash
./build/fem_solver examples/classic_bar_tension.inp output/classic_bar.vtk 1.0 dense
./build/fem_solver examples/classic_beam_tip_load.inp output/classic_beam.vtk 1.0 dense
./build/fem_solver examples/classic_solid_patch.inp output/classic_solid.vtk 1.0 dense
```

### 5.2 扩展关键字验证

```bash
./build/fem_solver --check-inp examples/advanced_keywords.inp
```

### 5.3 复杂功能案例（百级单元 + 接触）

```bash
./build/fem_solver examples/solid_block_200e.inp output/solid_200e.vtk 1.0 dense
./build/fem_solver examples/contact_two_blocks_144e.inp output/contact_144e.vtk 1.0 dense
```

说明：
- `solid_block_200e.inp`：200 个 C3D8 单元，可用于中等规模性能与稳定性测试。
- `contact_two_blocks_144e.inp`：上下块接触（含摩擦）非线性算例，会在终端打印增量/迭代残量过程。


### 5.4 大型案例（千级网格）

```bash
./build/fem_solver examples/cantilever_solid_2000e.inp output/cantilever_2000e.vtk 1.0 dense
```

### 5.5 FIELD OUTPUT 帧输出示例

```bash
./build/fem_solver examples/nonlinear_contact_field_output.inp output/nl_field.vtk 1.0 dense
```

当 `*Output, field, frequency=n` 或 `*Field Output, frequency=n` 生效时，会输出：
- 主结果：`output/nl_field.vtk`
- 帧结果：`output/nl_field_f0001.vtk`、`...`


---

## 6. 部署（可执行）

### 6.1 macOS（Apple Silicon / M4）

```bash
./scripts/deploy_macos_m4.sh
```

脚本会自动完成：平台检查、依赖安装（brew）、编译与冒烟测试。

### 6.2 Windows（PowerShell）

```powershell
.\scripts\deploy_windows.ps1
```

### 6.3 Windows（CMD）

```bat
scripts\deploy_windows.bat
```

---

## 7. 文档导航

- `docs/user_guide_cn.md`：用户使用帮助（命令、案例、FAQ）
- `docs/theory_manual_cn.md`：理论手册（离散方程、非线性、接触、约束、弧长）
- `docs/developer_manual_cn.md`：开发文档（架构、数据流、扩展点、调试建议）
- `docs/abaqus_static_alignment_plan.md`：后续对标 Abaqus 的分阶段路线图

---

## 8. 常见问题（简版）

1) **Singular matrix**
- 通常是约束不足（刚体自由度未被抑制）。

2) **非线性不收敛**
- 适当降低初始增量、增加 `*Controls` 迭代上限、调整接触罚系数与弧长参数。

3) **eigen/petsc 选了但不能用**
- 需要 CMake 检测到对应依赖；否则只能使用 dense/pcg。




## 19. Release 构建与 GitHub 发布

### 19.1 本地 Release 打包（当前平台）

```bash
./scripts/package_release_local.sh
```

该脚本会输出：`fem_solver-<os>-<arch>-release.tar.gz`。

### 19.2 GitHub 自动构建 macOS + Windows Release

仓库已提供工作流：`.github/workflows/release-build.yml`

触发方式：
- 推送 tag（如 `v0.4.0`）
- 或手动触发 `workflow_dispatch`

工作流会在 `macos-14` 与 `windows-2022` 上构建 Release 并上传产物：
- `fem_solver-macos-release.tar.gz`
- `fem_solver-windows-release.zip`

