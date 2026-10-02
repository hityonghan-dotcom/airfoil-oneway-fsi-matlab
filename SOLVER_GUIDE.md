# Solver guide / 求解器教学导读

This document is the teaching-facing interface for the MATLAB solver. It serves the same purpose as a C/C++ header file plus an algorithm note: students can learn the input, output, equations, and numerical sequence without first reading every implementation line.

本文档是面向教学的求解器接口说明。它相当于 C/C++ 的头文件加算法说明：学生无需先读完实现细节，也能知道输入、输出、控制方程与计算顺序。

## 1. Student entry points / 学生入口

| Goal / 目标 | File / 文件 | What students edit or run / 学生操作 |
| --- | --- | --- |
| One physical case / 单一物理工况 | **student_case_config.m** | Edit the labelled blocks A-E only / 只改 A-E 标注区 |
| Run one case / 运行单工况 | **run_student_case.m** | Run: r = run_student_case('demo') |
| Sweep incidence / 扫描攻角 | **lesson_02_angle_sweep.m** | Edit cfg.alphaDeg, then run the file |
| Public solver contract / 公开求解器接口 | **airfoil_solver.m** | Read the help text; normally do not edit |
| CFD implementation / CFD 实现 | **solve_airfoil_ns.m** | Read in the five numbered sections; advanced study only |
| One-way dynamics / 单向结构动力学 | **wing_dynamics.m** | Change stiffness, mass, or damping through configuration |

A beginner should only edit **student_case_config.m**. The run script creates the grid preview, solves the fixed-geometry flow, transfers loads to the structure, and saves the MATLAB result.

初学者只需修改 **student_case_config.m**。运行脚本会生成网格预览、计算固定翼型流场、把气动力传递给结构，并保存 MATLAB 结果。

### Minimal run / 最小运行方式

~~~matlab
% 1. Edit student_case_config.m.
% 2. Run one case.
r = run_student_case('demo');

% 3. Inspect physical outputs.
r.CL
r.CD
r.structure.tipBending
r.dynamics.theta
~~~

The default flow boundary condition is a uniform horizontal inlet: u = U_inf, v = 0. The top and bottom impose zero normal velocity, and the outlet extrapolates velocity with a reference pressure of zero. The course case uses a stationary airfoil at a prescribed geometric incidence.

默认流动边界条件为均匀水平来流：u = U_inf、v = 0。上下边界的法向速度为零，出口速度外推并取参考压力为零。课程案例采用固定翼型与给定几何攻角。

## 2. Public interface, like a header / 类似头文件的公开接口

### cfg input / 输入参数

| Group / 类别 | Important fields / 重要变量 | Meaning / 含义 |
| --- | --- | --- |
| Flow / 流动 | U, Re, nu, rho | Inlet speed, Reynolds number, viscosity, density / 来流速度、雷诺数、黏性、密度 |
| Geometry / 几何 | c, thickness, alphaDeg, center | Chord, NACA 00xx thickness, incidence, section position / 弦长、厚度、攻角、位置 |
| Grid / 网格 | Lx, Ly, nx, ny | Cartesian domain and pressure-cell count / 笛卡尔域与压力单元数 |
| Numerical controls / 数值控制 | dtScale, etaRatio, endTime, avgStart | Time step scale, penalty strength, transient duration, averaging start / 时间步、惩罚强度、总时长、平均开始时刻 |
| Structure / 结构 | span, EI, GJ, massPerLength, dampingRatio | Cantilever and reduced modal parameters / 悬臂梁与降阶模态参数 |

Keep dx = Lx/nx equal to dy = Ly/ny. If students change U, c, or Re while using the Reynolds-number definition, they must update nu = U*c/Re.

请保持 dx = Lx/nx 与 dy = Ly/ny 相等。如果学生改变 U、c 或 Re，并且仍采用雷诺数定义，则必须更新 nu = U*c/Re。

### r output / 输出结果

| Field / 变量 | Meaning / 含义 |
| --- | --- |
| history | Instantaneous time, CL, CD, CM, divergence, pressure residual, solid slip, CFL, speed, and lift per span / 瞬态载荷与数值诊断 |
| CL, CD, CM | Mean coefficients over [avgStart, endTime] / 平均时间窗内的气动系数 |
| U, V, p | Late-time mean velocity and pressure fields / 后期平均速度与压力场 |
| chiU, chiV | Brinkman solid fractions on staggered velocity faces / 交错速度面上的固体体积分数 |
| structure | Static cantilever response to mean CFD loads / 平均 CFD 载荷下的静态结构形变 |
| dynamics | Transient bending and torsion driven by the CFD history / CFD 瞬态载荷驱动的弯曲与扭转 |

## 3. Model scope / 模型范围

The package solves two-dimensional laminar incompressible flow around a fixed NACA 00xx section. It repeats the resulting two-dimensional load along a teaching cantilever span. It can illustrate low-Reynolds-number separation trends and a load-to-response calculation. It does not model turbulence, transition, finite-span physics, moving-grid feedback, or experimentally validated stall.

本程序求解固定 NACA 00xx 截面周围的二维层流不可压缩流动，并把二维截面载荷沿教学用悬臂翼展重复。它可用于展示低雷诺数分离趋势和“载荷到响应”的计算，不包含湍流、转捩、有限翼展效应、动网格反馈或实验验证的失速预测。

## 4. Governing equations / 控制方程

In the fluid region, the dimensional incompressible equations are:

在流体区域，程序使用以下有量纲不可压缩方程：

\[
\nabla\cdot\mathbf{u}=0,
\qquad
\frac{\partial\mathbf{u}}{\partial t}
+\nabla\cdot(\mathbf{u}\otimes\mathbf{u})
=-\frac{1}{\rho}\nabla p+\nu\nabla^2\mathbf{u}.
\]

The Cartesian-grid extension introduces a fixed Brinkman mask, chi. The body velocity is zero because the CFD geometry is fixed.

在笛卡尔网格上，程序加入固定的 Brinkman 掩膜 chi。由于 CFD 几何固定，刚体速度为零。

\[
\frac{\partial\mathbf{u}}{\partial t}
+\nabla\cdot(\mathbf{u}\otimes\mathbf{u})
=-\frac{1}{\rho}\nabla p+\nu\nabla^2\mathbf{u}
-\frac{\chi}{\eta}(\mathbf{u}-\mathbf{u}_s),
\qquad \mathbf{u}_s=0.
\]

Here chi = 0 in fluid and approaches 1 inside the airfoil. Eta is the penalty relaxation time. Smaller eta enforces no slip more strongly, but it may make the calculation less tolerant of numerical choices.

其中，流体区 chi = 0，翼型内部 chi 接近 1。eta 是惩罚松弛时间。较小的 eta 会更强地满足无滑移，但也会降低数值计算的宽容度。

## 5. Discretisation and projection method / 离散与投影法

### MAC grid / 交错网格

Pressure is stored at cell centres. Horizontal velocity u is stored on vertical faces, and vertical velocity v on horizontal faces. This staggered placement avoids the simplest pressure-velocity checkerboard pattern.

压力储存在单元中心；水平速度 u 储存在竖直单元面；竖直速度 v 储存在水平单元面。这种交错布置避免了最基本的压力与速度棋盘格问题。

### One time step / 一个时间步

The solver uses the following sequence in solve_airfoil_ns.m:

求解器在 solve_airfoil_ns.m 中按下列顺序运行：

1. **Momentum prediction / 动量预测**
   Explicit conservative convection and central-difference diffusion form u_hat, v_hat.

   显式保守形式对流项与中心差分黏性项得到 u_hat、v_hat。

2. **Implicit Brinkman step / 隐式 Brinkman 步**
   At each staggered face:
   \[
   u^*=\frac{u_{\mathrm{hat}}}{1+\Delta t\chi_u/\eta},
   \qquad
   v^*=\frac{v_{\mathrm{hat}}}{1+\Delta t\chi_v/\eta}.
   \]
   In fluid chi = 0 and this factor equals one. In the body, it damps velocity toward zero.

   在流体区 chi = 0，该因子为一；在固体内部，它把速度衰减至零。

3. **Pressure projection / 压力投影**
   The pressure equation uses the same Brinkman factor a:
   \[
   -\nabla\cdot(a\nabla p)
   =-\frac{\rho}{\Delta t}\nabla\cdot\mathbf{u}^*,
   \qquad
   a=\frac{1}{1+\Delta t\chi/\eta}.
   \]
   Then the staggered correction is:
   \[
   \mathbf{u}^{n+1}=\mathbf{u}^*
   -\frac{\Delta t}{\rho}a\nabla p.
   \]
   pressure_matrix assembles this sparse variable-coefficient operator and uses a Cholesky factorisation.

   然后在交错速度面上作压力修正。pressure_matrix 组装稀疏变系数矩阵，并采用 Cholesky 分解求解。

4. **Diagnostics / 数值诊断**
   The solver records dimensionless divergence, pressure residual, velocity inside nearly solid faces, and CFL number. Students should check these before interpreting force differences as physics.

   程序记录无量纲散度、压力残差、近乎实体面内的速度和 CFL 数。学生应先检查这些量，再讨论载荷变化是否反映物理。

The implementation uses first-order time advancement for the explicit momentum stage. It is intentionally compact for teaching rather than a production CFD method.

显式动量阶段采用一阶时间推进。它的目标是教学上的紧凑与可读性，而不是工业级 CFD 精度。

## 6. Brinkman immersed boundary / Brinkman 浸没边界法

airfoil_geometry.m creates the rotated NACA 00xx polygon. fraction in solve_airfoil_ns.m samples every velocity face with subcells by subcells points. It stores the fraction of a face covered by the airfoil in chiU or chiV.

airfoil_geometry.m 生成旋转后的 NACA 00xx 多边形。solve_airfoil_ns.m 中的 fraction 对每个速度面做 subcells x subcells 采样，并把翼型覆盖比例存入 chiU 或 chiV。

This is an immersed-boundary volume-penalisation method:

- The Cartesian mesh does not conform to the airfoil surface.
- The mask adds a resistance term inside the virtual solid.
- plot_grid_preview.m displays the mesh, airfoil, mask contour, and u/v face locations before marching.

这是一种浸没边界体积惩罚法：

- 笛卡尔网格不贴合翼型表面。
- 掩膜在虚拟固体内部增加阻尼项。
- plot_grid_preview.m 在时间推进前显示网格、翼型、掩膜等值线和 u/v 面位置。

## 7. Force extraction and one-way FSI / 载荷提取与单向流固耦合

The penalty force acting on the virtual fluid is:

\[
\mathbf{f}_{penalty}
=-\frac{\rho}{\eta}\chi\mathbf{u}^{n+1}.
\]

The body receives the opposite reaction. The code also adds the virtual-fluid momentum correction before integrating force and moment:

\[
D'=\sum q_x\Delta A,\qquad
L'=\sum q_y\Delta A,\qquad
M'=\sum\left[(y-y_e)q_x-(x-x_e)q_y\right]\Delta A.
\]

wing_dynamics.m reads the already-computed L'(t) and M'(t) and solves:

\[
M\ddot{\mathbf q}+C\dot{\mathbf q}+K\mathbf q=\mathbf Q(t),
\qquad \mathbf q=(h,\theta)^T.
\]

newmark_linear.m uses average-acceleration Newmark integration with beta = 1/4 and gamma = 1/2.

wing_dynamics.m 读取已经算好的 L'(t) 与 M'(t)，求解上式。newmark_linear.m 使用平均加速度 Newmark 积分，参数为 beta = 1/4、gamma = 1/2。

### Coupling boundary / 耦合边界

\[
\text{fixed CFD geometry}
\longrightarrow
\{L'(t),D'(t),M'(t)\}
\longrightarrow
\{h(t),\theta(t)\}.
\]

No h or theta enters airfoil_geometry, the mask, or the CFD boundary conditions. The moving airfoil in the motion video is a deliberately magnified visual overlay, not a moving CFD boundary.

h 和 theta 不会进入 airfoil_geometry、掩膜或 CFD 边界条件。运动视频中的翼型是经过刻意放大的可视化叠加，并不是运动 CFD 边界。

## 8. Suggested student investigations / 建议学生开展的物理问题

| Question / 问题 | Change / 修改 | Compare / 比较 |
| --- | --- | --- |
| Incidence trend / 攻角趋势 | alphaDeg | Mean CL, CD, wake low-speed region |
| Reynolds-number trend / 雷诺数趋势 | Set Re, then update nu | Load history and separation trend |
| Grid sensitivity / 网格敏感性 | Scale nx, ny while keeping dx=dy | CL, CD, maxDiv, windowChange |
| Penalty sensitivity / 惩罚参数敏感性 | etaRatio | Solid slip and forces |
| Structural stiffness / 结构刚度 | EI, GJ | tipBending, tipTwistDeg, dynamics |
| Damping and resonance / 阻尼与响应 | dampingRatio, mass fields | Oscillation amplitude and settling |

Students should label all aerodynamic conclusions as numerical low-Reynolds-number trends. The package does not establish a validated physical stall angle.

学生应把所有气动结论表述为低雷诺数数值趋势。本程序不能给出经过验证的真实失速攻角。

## 9. Lecture slide plan / 教学 PPT 建议

| Slide / 页 | Chinese teaching point / 中文讲授重点 | English title and content cue / 英文标题与内容提示 |
| --- | --- | --- |
| 1 | 课程问题：固定翼型流场如何产生结构响应 | **Airfoil flow and one-way structural response**. Show the end-to-end workflow. |
| 2 | 明确模型范围：二维、层流、低 Re、固定几何 | **Model scope**. State what is included and omitted. |
| 3 | 让学生先运行单工况并看网格预览 | **Student workflow**. Show student_case_config.m and one run command. |
| 4 | 介绍 MAC 交错网格与变量位置 | **Staggered Cartesian grid**. Use the generated grid preview. |
| 5 | 写出连续 NS 方程和边界条件 | **Incompressible Navier-Stokes equations**. Define u, p, rho, nu. |
| 6 | 逐步展示预测、惩罚、投影 | **One CFD time step**. Show the three equations in Section 5. |
| 7 | 解释 chi、eta 和浸没边界 | **Brinkman immersed boundary**. Use the mask zoom from the grid preview. |
| 8 | 由惩罚反力得到升力、阻力和力矩 | **Loads from the penalty reaction**. Connect qx, qy to L', D', M'. |
| 9 | 明确单向耦合箭头及 Newmark 结构响应 | **One-way load transfer**. Show that motion does not enter CFD. |
| 10 | 让学生做参数研究并讨论数值可信度 | **Student investigations and numerical checks**. Use CL, CD, divergence, slip, and sensitivity results. |

For slides 4, 7, and 9, use files generated by the package rather than illustrative stock images: grid_preview_alpha_*.png, flow_field.png, and oneway_motion_alpha_*.mp4.

第 4、7、9 页建议直接使用程序生成的 grid_preview_alpha_*.png、flow_field.png 与 oneway_motion_alpha_*.mp4，这样图、方程和代码能一一对应。

## 10. Reading order for advanced students / 进阶学生的阅读顺序

1. student_case_config.m
2. airfoil_solver.m
3. solve_airfoil_ns.m, only its five numbered blocks
4. plot_grid_preview.m
5. wing_structure.m, wing_dynamics.m, then newmark_linear.m
6. run_sensitivity.m

This order keeps the physical model visible before students encounter implementation detail.

该顺序先让学生建立物理模型，再接触具体实现细节。
