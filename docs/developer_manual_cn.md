# mini_fem_solver 开发者教材（架构 + 算法 + 实施）

> 目标：让开发者在不了解历史上下文的情况下，也能独立完成“新增功能—验证—交付”。

---

## 1. 开发目标与原则

### 1.1 目标
- 保持可运行。
- 保持可读。
- 保持可扩展。
- 保持跨平台。

### 1.2 原则
- 先正确后优化。
- 先最小可验证再扩展。
- 每次改动必须有可复现验证命令。

---

## 2. 架构总览

### 2.1 模块划分
- `main.cpp`：CLI 与任务入口。
- `inp_parser.cpp`：输入解析。
- `solver.cpp`：核心算法。
- `vtk_writer.cpp`：结果导出。

### 2.2 数据流
1. CLI 读取参数。
2. parser 生成 `Model`。
3. solver 计算 `Result`。
4. writer 输出 VTK。

---

## 3. 核心数据结构（必须读懂）

### 3.1 `Model`
承载网格、材料、载荷、约束、接触、分析步参数。

### 3.2 `Step`
包含：
- 迭代上限；
- 容差；
- cutback 参数；
- 弧长参数。

### 3.3 `Result`
位移、反力、单元结果字段。

### 3.4 `MpcConstraint`
显式多点约束结构。

---

## 4. parser 实现细节（`src/inp_parser.cpp`）

### 4.1 状态机实现
通过 `currentSection` 区分当前关键字上下文。

### 4.2 关键字解析流程
1. `parseKeyword` 分离主关键字与参数。
2. 根据关键字更新 parser 状态。
3. 读数据行填充 `Model`。
4. 文件结束时做一致性校验。

### 4.3 已扩展关键字
- `*MPC`
- 扩展 `*CONTROLS`
- `*Contact Pair` 参数化读取

### 4.4 parser 扩展模板（实操）
1. 加关键字识别分支。
2. 加数据行解析分支。
3. 加字段校验。
4. 加示例输入。
5. 更新文档和测试命令。

---

## 5. solver 总体框架（`src/solver.cpp`）

### 5.1 函数职责分层
- 元素级：`addTruss2/addBeam31/addShell4/addC3D8`
- 组装级：`assemble`
- 约束/接触：`applyBC/applyCouplingPenalty/applyContactPenalty`
- 线性求解：`solveLinearSystem`
- 非线性流程：`solveNonlinear`

### 5.2 线性求解路径
- 组装 K、f
- 施加边界
- 解线性系统
- 计算反力与单元结果

### 5.3 非线性求解路径
- 外层载荷增量循环
- 内层 Newton 迭代循环
- 失败触发 cutback
- 可选弧长参数调节

---

## 6. 关键算法详解

### 6.1 节点索引缓存（性能点）
`nodeIdx` 采用缓存映射，避免频繁 O(N) 节点查找。

### 6.2 外载统一构造
`buildExternalLoad` 统一处理节点载荷 + DLOAD/DSLOAD，减少重复代码，提升可维护性。

### 6.3 约束增广系统
`solveWithMpcLagrange` 将 Coupling 与 MPC 统一为线性约束行，组装增广系统解增量。

### 6.4 接触与摩擦
- 法向：穿透时激活罚力和法向切线。
- 切向：按库仑阈值决定切向刚度投影。

### 6.5 弧长+cutback
- 收敛快时增大步长；
- 发散时减小步长并回退；
- 弧长半径有上下限。

---

## 7. 可读性与可扩展性实践建议

### 7.1 命名
- 用“动作 + 对象”命名函数。
- 用“物理量 + 单位语义”命名变量。

### 7.2 注释
写“为什么这样做”，而不仅是“做了什么”。

### 7.3 解耦
- parser 只负责“读 + 校验”，不做求解逻辑。
- solver 不依赖 CLI 细节。
- writer 不依赖 solver 内部临时变量。

### 7.4 扩展路径
新增算法优先以“独立函数 + 最小调用点”的方式落地。

---

## 8. 性能优化建议（按优先级）

1. 降低热点 O(N²) 查找。
2. 减少重复装配与重复解析。
3. 引入稀疏矩阵存储（后续重点）。
4. 对大循环做 OpenMP 并行。
5. 在保持正确性的前提下做内存局部性优化。

---

## 9. 稳定性与调试手册

### 9.1 常见报错与诊断
- `Singular matrix`：约束不足或冲突。
- `Newton did not converge`：步长/参数过激。
- 接触振荡：罚系数过大、摩擦参数不当。

### 9.2 推荐排查顺序
1. 先跑最小案例。
2. 检查边界与载荷单位。
3. 关闭摩擦仅保留法向接触。
4. 降低载荷增量与提高迭代次数。

---

## 10. 测试与回归策略

### 10.1 基础构建
```bash
cmake -S . -B build
cmake --build build -j
```

### 10.2 CLI 功能验证
```bash
./build/fem_solver --help
./build/fem_solver --solvers
./build/fem_solver --check-inp examples/advanced_keywords.inp
```

### 10.3 经典回归
```bash
./build/fem_solver examples/classic_bar_tension.inp output/classic_bar.vtk 1.0 dense
./build/fem_solver examples/classic_beam_tip_load.inp output/classic_beam.vtk 1.0 dense
./build/fem_solver examples/classic_solid_patch.inp output/classic_solid.vtk 1.0 dense
```

---

## 11. 新功能开发流程（可直接照做）

1. 先写一份最小输入案例。
2. 在 parser 增加关键字支持。
3. 在模型结构中增加字段。
4. 在 solver 增加装配和求解逻辑。
5. 跑最小案例验证。
6. 回归已有案例。
7. 更新 README + 文档。
8. 提交 PR。

---

## 12. 算法级伪代码模板

### 12.1 非线性主循环
```text
for each load increment:
    build target load
    for each Newton iteration:
        assemble K and fint
        R = fext - fint
        if converged: break
        solve K * du = R (with constraints)
        u += du
    if not converged:
        cutback and retry
```

### 12.2 接触装配
```text
for each contact pair:
    compute normal n and gap
    if penetration:
        add normal penalty tangent
        add normal residual
        compute tangential slip
        apply Coulomb limit
        add tangential tangent/residual
```

### 12.3 MPC/Lagrange
```text
collect linear constraints C*u = g
assemble augmented matrix [K C^T; C 0]
solve for [du, dlambda]
```

---

## 13. 代码审查清单（Reviewer Checklist）

- 是否引入了新物理假设？是否文档化？
- 是否提供最小案例？
- 是否有回归命令？
- 是否破坏了现有 CLI 行为？
- 是否添加了必要中文注释？
- 是否考虑跨平台编译？

---

## 14. 长期路线（面向 Abaqus 对标）

### 14.1 数值层
- 更标准弧长联立算法。
- 稀疏矩阵 + 预条件。
- 更完整接触算法（surface-to-surface）。

### 14.2 物理层
- 更完整塑性/损伤模型。
- 壳与实体的多积分点统一管理。
- 更系统的 hourglass 控制。

### 14.3 工程层
- 自动回归脚本。
- 基准算例数据库。
- 文档与代码双向追踪。

---

## 附录 A：新成员两周上手计划

- 第 1~2 天：读 README + 跑经典案例。
- 第 3~5 天：读 parser + 修改一个关键字。
- 第 6~8 天：读 solver + 加一个小输出字段。
- 第 9~11 天：实现一个小算法改进。
- 第 12~14 天：完成回归 + 文档 + PR。

---

## 附录 B：常见“踩坑”速查

- 单位混乱（N、mm、MPa）
- 约束重复或冲突
- 罚系数数量级过大
- 步长过激导致迭代发散
- 把“看起来收敛”误认为“物理正确”

