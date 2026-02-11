# mini_fem_solver 开发文档（中文）

> 本文面向二次开发者，强调代码结构、数据流、扩展路径与调试建议。

---

## 1. 目录结构

- `include/`
  - `fem_types.hpp`：核心数据结构（Model/Element/Step/Result 等）
  - `inp_parser.hpp` / `solver.hpp` / `vtk_writer.hpp`
- `src/`
  - `inp_parser.cpp`：`.inp` 关键字解析与模型构建
  - `solver.cpp`：组装、线性/非线性求解、接触与约束
  - `vtk_writer.cpp`：结果输出
  - `main.cpp`：CLI 交互入口
- `examples/`：回归与演示输入
- `scripts/`：部署脚本
- `docs/`：用户/理论/路线图文档

---

## 2. 核心数据模型

建议先读 `include/fem_types.hpp`：

- `Model`：全局容器
  - 网格（nodes/elements）
  - 集合（nsets/elsets）
  - 材料、边界、载荷
  - 接触、耦合、MPC
  - 分析步参数 `Step`
- `Step`：求解控制
  - 非线性迭代、容差、cutback
  - 弧长半径及调节参数
- `Result`：位移/反力/单元场输出

---

## 3. 程序主数据流

1. `main.cpp`
   - 解析 CLI
   - 读取 `.inp` -> `Model`
   - 调用 `solve(model)`
   - 写出 VTK
2. `inp_parser.cpp`
   - 逐行扫描关键字
   - 填充 `Model` 字段
   - 末尾执行一致性校验（节点引用、MPC 引用）
3. `solver.cpp`
   - `assemble()` 完成全局 K 与内力
   - 线性/非线性分支
   - 非线性中反复组装、线性化、更新

---

## 4. 扩展关键字的标准步骤

以新增 `*XXX` 为例：

1. 在 `inp_parser.cpp` 的 keyword 分支加入识别。
2. 在数据行分支读取数值并写入 `Model`。
3. 若涉及求解，在 `solver.cpp` 中增加装配逻辑。
4. 增加最小算例到 `examples/`。
5. 在 `README.md` 与用户文档补充说明。

---

## 5. 扩展一个新单元类型的建议流程

1. 在 `ElementType` 枚举加入类型。
2. 在 parser 的 `parseElementType()` 增加映射。
3. 在 solver 中新增 `addXXX()`：
   - 构造 DOF map
   - 计算 `ke`、`fint`
   - `addToGlobal()` 组装
4. 增加回归算例验证刚度与位移方向。

---

## 6. 非线性收敛调试建议

- 打印每步残量范数与增量范数（可临时加日志）。
- 从最小模型开始（1~2 单元）定位问题。
- 先关摩擦再开摩擦，分离接触问题来源。
- 调低接触罚系数与初始增量，观察收敛敏感性。
- 关注边界约束是否充分，避免奇异刚度矩阵。

---

## 7. 跨平台构建建议

- CMake 选项：
  - `FEM_ENABLE_EIGEN`
  - `FEM_ENABLE_PETSC`
  - `FEM_ENABLE_OPENMP`
- 优先保证 `dense` 后端始终可运行，作为兜底路径。
- Windows 推荐 Ninja + VS 工具链。

---

## 8. 回归检查清单（建议 PR 前执行）

```bash
cmake -S . -B build
cmake --build build -j
./build/fem_solver --help
./build/fem_solver --examples
./build/fem_solver --check-inp examples/advanced_keywords.inp
./build/fem_solver examples/classic_bar_tension.inp output/classic_bar.vtk 1.0 dense
./build/fem_solver examples/classic_beam_tip_load.inp output/classic_beam.vtk 1.0 dense
./build/fem_solver examples/classic_solid_patch.inp output/classic_solid.vtk 1.0 dense
```

---

## 9. Roadmap 对接

开发优先级可参照：
- `docs/abaqus_static_alignment_plan.md`

建议采用“功能分支 + 最小案例 + 文档同步 + 回归脚本”节奏迭代。

