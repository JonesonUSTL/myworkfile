# mini_fem_solver 理论手册（中文）

> 本文面向“会用有限元但想看清实现骨架”的读者，重点解释当前程序中的离散方程、非线性迭代、接触与约束处理逻辑。

---

## 1. 总体离散方程

静力平衡方程离散后可写为：

\[
\mathbf{R}(\mathbf{u},\lambda)=\mathbf{f}_{ext}(\lambda)-\mathbf{f}_{int}(\mathbf{u})=\mathbf{0}
\]

其中：
- \(\mathbf{u}\)：全局自由度向量
- \(\lambda\)：载荷比例参数
- \(\mathbf{f}_{ext}\)：外载向量（含幅值与体载近似转换）
- \(\mathbf{f}_{int}\)：内力向量（由单元应力积分得到）

牛顿迭代线性化：

\[
\mathbf{K}_t \Delta\mathbf{u} = \mathbf{R}
\]

\(\mathbf{K}_t\) 为切线刚度（材料刚度 + 几何/接触/约束贡献）。

---

## 2. 单元层与组装

程序统一采用“单元局部矩阵/向量 -> 全局组装”流程：

1. 读取单元节点坐标、材料参数、状态变量。
2. 在积分点计算应变、应力与切线刚度（J2 return-mapping 用于实体）。
3. 积分得到单元刚度 \(\mathbf{k}_e\) 与内力 \(\mathbf{f}_{int,e}\)。
4. 按 DOF 映射装配到全局 \(\mathbf{K}\)、\(\mathbf{f}_{int}\)。

---

## 3. 材料本构简述

### 3.1 线弹性

采用各向同性线弹性矩阵 \(\mathbf{D}(E,\nu)\)：

\[
\boldsymbol{\sigma}=\mathbf{D}\,\boldsymbol{\varepsilon}
\]

### 3.2 双线性弹塑性（杆）

通过屈服应变阈值判断弹性/塑性区，塑性段切线模量使用硬化模量 \(H\)。

### 3.3 J2 塑性（实体）

采用试应力 + 回归映射（return-mapping）更新应力与一致切线近似。

---

## 4. 约束方程（Coupling + MPC）

程序将约束统一为：

\[
\sum_i c_i u_i = g
\]

并采用 Lagrange 乘子增广系统：

\[
\begin{bmatrix}
\mathbf{K} & \mathbf{C}^T \\
\mathbf{C} & \mathbf{0}
\end{bmatrix}
\begin{bmatrix}
\Delta\mathbf{u} \\
\Delta\boldsymbol{\lambda}_c
\end{bmatrix}
=
-\begin{bmatrix}
\mathbf{R} \\
\mathbf{g}
\end{bmatrix}
\]

优点：
- 约束表达直观；
- 易于扩展多主节点 MPC、后续接触拉格朗日框架。

---

## 5. 接触与摩擦（罚函数）

当前采用简化 node-to-node 配对近似：

1. 计算主从节点相对方向法向 \(\mathbf{n}\)。
2. 法向间隙 \(g_n = (u_s-u_m)\cdot n\)。
3. 当 \(g_n<0\)（压入）时激活法向罚力：

\[
p_n = k_n g_n
\]

切线贡献：

\[
\mathbf{K}_n \sim k_n (\mathbf{n}\otimes\mathbf{n})
\]

摩擦切向：
- 切向相对位移 \(\mathbf{u}_t\)
- 试算切向力 \(t_{trial}=k_t\|u_t\|\)
- 若 \(t_{trial}>\mu f_n\) 则按库仑上限退化切向刚度。

---

## 6. 非线性步进、cutback 与弧长参数

### 6.1 常规增量 Newton

每个载荷增量内迭代到残量范数收敛：
- 收敛快：可放大下步增量
- 收敛差/失败：减小增量（cutback）

### 6.2 弧长参数化

程序当前实现“工程可用的弧长半径调节框架”：
- `arcLengthRadius`
- `arcLengthGrowFactor` / `arcLengthShrinkFactor`
- `arcLengthMinRadius` / `arcLengthMaxRadius`

其作用是控制每步增量尺度，提升极限点附近路径追踪稳定性。

---

## 7. 线性求解器后端

- `dense`：致密高斯消元，依赖最少。
- `eigen`：Eigen 稀疏因子分解（需依赖）。
- `petsc`：PETSc KSP（需依赖）。
- `pcg`：OpenMP 并行 CG（适用于近 SPD 问题）。

---

## 8. 结果输出（VTK）

程序输出：
- 节点位移
- 节点反力
- 元素轴力/应变/应力/von Mises（按单元类型可用性）

推荐用 ParaView 进行：
- Warp By Vector（查看变形）
- 着色查看应力应变场

---

## 9. 当前局限与下一步

- 接触为简化配对模型，尚非完整 surface-to-surface。
- 壳与实体多积分点、严格 hourglass 控制仍可继续增强。
- 弧长法仍为工程化实现，后续可扩展更标准 Riks 联立方程。

