# mini_fem_solver 使用帮助（中文）

## 1. 工具定位

`mini_fem_solver` 是一个可在 macOS/Windows/Linux 构建运行的 Abaqus 风格 `inp` 静力学求解器，支持：
- 线性/非线性静力（Newton-Raphson）
- Truss / Beam / Shell / Solid 基础单元
- 接触罚函数（法向 + 摩擦）
- Coupling/Kinematic + MPC 显式约束方程（Lagrange 增广）
- VTK 结果输出（位移/反力/应力相关字段）

## 2. 命令行交互

```bash
./build/fem_solver --help
./build/fem_solver --examples
./build/fem_solver <input.inp> <output.vtk> [deform_scale] [solver=dense|eigen|petsc|pcg]
```

参数说明：
- `input.inp`：Abaqus 风格输入文件。
- `output.vtk`：输出路径，目录会自动创建。
- `deform_scale`：变形显示缩放，默认 `1.0`。
- `solver`：线性求解后端，默认 `dense`。

## 3. 推荐快速验证案例

```bash
./build/fem_solver examples/classic_bar_tension.inp output/classic_bar.vtk 1.0 dense
./build/fem_solver examples/classic_beam_tip_load.inp output/classic_beam.vtk 1.0 dense
./build/fem_solver examples/classic_solid_patch.inp output/classic_solid.vtk 1.0 dense
```

## 4. 关键字提示

当前重点支持：
- 网格与集合：`*Node` `*Element` `*Nset` `*Elset`
- 材料与截面：`*Material` `*Elastic` `*Plastic` `*Truss Section` `*Beam Section` `*Shell Section`
- 边界与载荷：`*Boundary` `*Cload` `*Dload` `*Dsload` `*Amplitude`
- 非线性扩展：`*Step` `*Static, riks` `*Controls`
- 约束/接触：`*Coupling` `*Kinematic` `*MPC` `*Contact Pair`

## 5. 常见问题

1) **报错 Singular matrix**
- 检查是否约束不足（刚体位移/转动未抑制）。

2) **非线性不收敛**
- 适当减小初始增量，增加 `*Controls` 的迭代上限，或调整接触罚系数。

3) **PETSc/Eigen 后端无法使用**
- 确认系统已安装对应依赖，并在 CMake 配置时开启选项。

## 6. 一键部署

- macOS（M4/ARM）：`scripts/deploy_macos_m4.sh`
- Windows（PowerShell）：`scripts/deploy_windows.ps1`
- Windows（CMD）：`scripts/deploy_windows.bat`
