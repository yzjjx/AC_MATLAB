# Y_a_reference 与 Y_a 的雅可比计算流程

## 1. 新流程与运行顺序

当前 [generate_Y_a.m](../MATLAB/generate_Y_a.m) 直接使用前置脚本提取的矩阵，完成“构造参考力矩、更换参数、对独立参数求雅可比”的流程。不再依赖原始回归矩阵 `Y` 或 `generate_Y_pi.m`。

用户指定的完整原始模型运行顺序：

```matlab
addpath('MATLAB');
RNEA_SYM;
Mat_out;
add_friction;
generate_Y_a;
```

沿用原平面二连杆七参数模型时：

```matlab
addpath('MATLAB');
RNEA_SYM;
data_verify;
Mat_out;
add_friction;
generate_Y_a;
```

`Mat_out` 必须在 `add_friction` 之前运行，因为它提取的是无摩擦刚体动力学矩阵。

完整原始 `RNEA_SYM` 保留质心 y 方向分量，当前二连杆生成 10 个参数组合；`data_verify` 将质心设置在连杆 x 轴后，生成原七参数组合。原硬编码七列的 `AC_CON` 仅适用于对应的七参数模型。

## 2. 前置脚本如何提取 M、C、G

[Mat_out.m](../MATLAB/Mat_out.m) 使用无摩擦逆动力学：

$$
\tau_{\mathrm{rigid}}=M(q)\ddot q+C(q,\dot q)\dot q+G(q).
$$

惯性矩阵是力矩对实际加速度的雅可比：

$$
M(q)=\frac{\partial\tau_{\mathrm{rigid}}}{\partial\ddot q}.
$$

```matlab
M = simplify(jacobian(tau_input,mat_ddq));
```

令速度和加速度为零，得到当前模型的重力项：

$$
G(q)=\tau_{\mathrm{rigid}}(q,0,0).
$$

```matlab
G = simplify(subs(tau_input, ...
    [mat_dq;mat_ddq],sym(zeros(2*N,1))));
```

科氏/离心矩阵使用 Christoffel 构造：

$$
C_{i\ell}(q,\dot q)=\frac12\sum_{k=1}^{N}
\left(
\frac{\partial M_{i\ell}}{\partial q_k}
+\frac{\partial M_{ik}}{\partial q_\ell}
-\frac{\partial M_{\ell k}}{\partial q_i}
\right)\dot q_k.
$$

`C` 中始终使用实际速度 `dq`。只有其右侧所乘的速度向量在参考补偿力矩中改为 `vr`。

## 3. 读取摩擦并核对前置模型

`add_friction` 保存原力矩为 `tau_before_friction`，然后增加摩擦：

$$
\tau=\tau_{\mathrm{before\_friction}}+F(q,\dot q).
$$

新脚本首先验证：

$$
\tau_{\mathrm{before\_friction}}-M\ddot q-C\dot q-G=0.
$$

不满足该式时，说明矩阵与本次力矩不一致，不能继续参数化。然后读取新增摩擦：

```matlab
F = simplify(expand(tau-tau_before_friction));
```

当前黏性摩擦为：

$$
F(\dot q)=\operatorname{diag}(f_{v1},\ldots,f_{vN})\dot q.
$$

一般附加速度项可以依赖位置，记为 `F(q,dq)`。脚本检查其不含实际加速度，并在位置保持不变时满足：

$$
F(q,-\dot q)=-F(q,\dot q).
$$

奇对称性是结构检查，不能独立证明任意残差都是物理摩擦。当前动力学采用 `+F` 补偿力矩约定，耗散摩擦还应满足：

$$
\dot q^TF(q,\dot q)\ge0.
$$

## 4. 构造参考量及原物理参数形式的力矩

跟踪误差为：

$$
e=q-q_d,\qquad \dot e=\dot q-\dot q_d.
$$

使用常矩阵增益：

$$
v_r=\dot q_d-\Lambda e,\qquad
a_r=\ddot q_d-\Lambda\dot e=\dot v_r,\qquad
s=\dot q-v_r=\dot e+\Lambda e.
$$

当前默认：

$$
\Lambda=3I_N.
$$

若增益随时间变化，参考速度的导数还包含额外项：

$$
\dot v_r=\ddot q_d-\Lambda\dot e-\dot\Lambda e.
$$

新脚本保留独立参考符号 `mat_vr、mat_ar`，并构造两种物理参数形式的补偿力矩：

```matlab
tau_a_reference_physical = simplify(M*mat_ar+C*mat_vr+G+F);
tau_a_physical = simplify(M*a_r+C*v_r+G+F);
```

前者的参考量是独立输入，后者已代入误差与期望轨迹定义。摩擦均保持实际 `dq`。

## 5. 为什么不能直接对质量、质心变量求雅可比

动力学中存在物理参数乘积，例如：

$$
m_1l_{c1}^2,\qquad m_2L_1l_{c2},\qquad m_2l_{c2}.
$$

直接对原始质量、质心等变量求导，通常不能得到满足齐次线性关系的普通 AC 回归矩阵。例如：

$$
\tau_a=ml_c^2a_r,
\qquad
\frac{\partial\tau_a}{\partial[m,l_c]}
=\begin{bmatrix}l_c^2a_r&2ml_ca_r\end{bmatrix}.
$$

上述雅可比再乘原变量向量，得到的不是原力矩：

$$
\begin{bmatrix}l_c^2a_r&2ml_ca_r\end{bmatrix}
\begin{bmatrix}m\\l_c\end{bmatrix}
=3ml_c^2a_r.
$$

应先定义固定参数组合，再使用独立符号表示：

$$
\pi_1=ml_c^2,\qquad
\tau_a=\pi_{a1}a_r,\qquad
\frac{\partial\tau_a}{\partial\pi_{a1}}=a_r.
$$

因此代码区分：

| 变量 | 含义 |
|---|---|
| `pi_param` | 固定参数组合的物理定义，包含质量、质心、惯量等符号 |
| `pi_a` | 与这些组合一一对应的独立求导符号 |
| `pi_hat` | 控制器中的数值参数估计，顺序须与本次 `pi_param` 相同 |

代码使用 `pi_a` 而非变量名 `pi`，避免覆盖 MATLAB 圆周率常数。

## 6. 自动更换参数的具体实现

脚本将参考补偿力矩和所有矩阵堆叠，确保新参数不仅能恢复力矩，还能恢复各动力学矩阵：

```matlab
physical_stack_a = [tau_a_reference_physical;M(:);C(:);G;F];
```

运动状态、参考量和默认已知的重力 `g` 不参与未知物理参数组合。辅助函数 `group_parameters_a` 用 `coeffs` 提取剩余物理符号的单项式与其系数向量。

例如，当前平面模型中以下三项在全部堆叠表达式中具有相同系数：

$$
I_{zz1},\qquad m_1l_{c1}^2,\qquad m_2L_1^2.
$$

它们被合并为：

$$
A=I_{zz1}+m_1l_{c1}^2+m_2L_1^2.
$$

若两项的系数向量成数值常数倍，也可合并：

$$
b_2=c\,b_1
\quad\Longrightarrow\quad
b_1\theta_1+b_2\theta_2=b_1(\theta_1+c\theta_2).
$$

这里的比例 `c` 不能依赖运动状态、参考量或物理符号。

该过程不使用随机采样、SVD/QR，仅合并相同或常数倍的系数。对一般模型，不保证已消除所有更复杂的线性依赖。因此，`pi_a` 是独立求导变量，不表示所有回归列都已证明独立、可辨识。

当前实现要求待估计物理参数保持符号形式，且改写后无与未知参数无关的非零偏置项；出现这种偏置时会报错。已知量配置在 `known_symbols_a` 中设置。

独立参数符号由以下代码生成：

```matlab
pi_a = sym('pi_a_',[parameter_count_a,1],'real');
```

对于原平面二连杆，如果七个物理定义完全匹配，脚本只调整参数排列顺序以兼容 `AC_CON`：

$$
\pi_{\mathrm{param}}=
\begin{bmatrix}f_{v1}\\A\\B\\H_1\\D\\H_2\\f_{v2}\end{bmatrix}.
$$

$$
A=I_{zz1}+m_1l_{c1}^2+m_2L_1^2,\qquad
B=m_2L_1l_{c2},\qquad H_1=m_1l_{c1}+m_2L_1.
$$

$$
D=I_{zz2}+m_2l_{c2}^2,\qquad H_2=m_2l_{c2}.
$$

更一般的模型按本次输出的物理定义与列序使用。

## 7. 构造新参数形式的 tau_a，直接求雅可比

将改写后的堆叠表达式拆回矩阵，得到依赖独立参数的 `M_a、C_a、G_a、F_a`。原物理形式 `M、C、G` 保留在工作区。

核心计算只有以下几步：

```matlab
tau_a_reference = simplify(M_a*mat_ar+C_a*mat_vr+G_a+F_a);
tau_a = simplify(M_a*a_r+C_a*v_r+G_a+F_a);

Y_a_reference = simplify(jacobian(tau_a_reference,pi_a));
Y_a = simplify(jacobian(tau_a,pi_a));
```

因为改写后的力矩对独立参数齐次线性：

$$
\tau_a=\sum_{j=1}^{p}y_{aj}\pi_{aj},
\qquad
\frac{\partial\tau_a}{\partial\pi_{aj}}=y_{aj}.
$$

所以：

$$
Y_a=\begin{bmatrix}y_{a1}&\cdots&y_{ap}\end{bmatrix},
\qquad
\tau_a=Y_a\pi_a.
$$

这直接实现了所要求的雅可比公式。

## 8. 原平面七参数模型的结果

记：

$$
c_2=\cos q_2,\qquad s_2=\sin q_2.
$$

新的求解流程得到的独立参考矩阵仍为：

$$
Y_{a,\mathrm{reference}}=
\begin{bmatrix}
\dot q_1&a_{r1}&\phi_{13}&g\sin q_1&a_{r1}+a_{r2}&g\sin(q_1+q_2)&0\\
0&0&\phi_{23}&0&a_{r1}+a_{r2}&g\sin(q_1+q_2)&\dot q_2
\end{bmatrix}.
$$

$$
\begin{aligned}
\phi_{13}&=2c_2a_{r1}+c_2a_{r2}
-\dot q_2s_2v_{r1}-(\dot q_1+\dot q_2)s_2v_{r2},\\
\phi_{23}&=c_2a_{r1}+\dot q_1s_2v_{r1}.
\end{aligned}
$$

第 3 列就是对独立耦合参数求导：

$$
\frac{\partial\tau_{a,\mathrm{reference}}}{\partial\pi_{a3}}
=\begin{bmatrix}\phi_{13}\\\phi_{23}\end{bmatrix}.
$$

摩擦列保留实际速度，科氏矩阵内部的速度也保留实际速度，不能对原表达式整体执行 `dq -> vr`。

## 9. 符号验证

脚本检查独立参数的线性回归恒等式：

$$
\tau_{a,\mathrm{reference}}-Y_{a,\mathrm{reference}}\pi_a=0,
\qquad \tau_a-Y_a\pi_a=0.
$$

代回物理参数定义，应恢复原物理形式：

$$
\left.\tau_a\right|_{\pi_a=\pi_{\mathrm{param}}}
-\tau_{a,\mathrm{physical}}=0.
$$

参考量取实际运动量时，应恢复完整原力矩：

$$
Y_{a,\mathrm{reference}}(q,\dot q,\dot q,\ddot q)
\pi_{\mathrm{param}}-\tau=0.
$$

参考矩阵代入误差定义后，应与直接对 `tau_a` 求雅可比的结果一致：

$$
\left.Y_{a,\mathrm{reference}}\right|_{v_r=\dot q_d-\Lambda e,\,
a_r=\ddot q_d-\Lambda\dot e}-Y_a=0.
$$

还检查 `Y_a` 不含 `pi_a`、未消除的未知物理符号或实际加速度 `ddq`。这些检查证明同一模型内部的代数一致性，不证明前置物理模型正确，也不证明在线参数估计满足持续激励或已经收敛。

## 10. 数值导出与矩阵还原

代入已知重力后可导出数值函数：

```matlab
Y_a_known = subs(Y_a,sym('g'),sym('9.81'));
Ya_fun = matlabFunction(Y_a_known, ...
    'Vars',{mat_q,mat_dq,mat_qd,mat_dqd,mat_ddqd});
```

控制器使用：

$$
\tau_{\mathrm{cmd}}=Y_a\hat\pi-K_s s.
$$

参数估计按本次 `pi_param` 的定义与列顺序初始化。

新 `M_a` 等矩阵使用 `pi_a`。若要继续用物理质量、质心、惯量代入数值，需要先恢复物理参数形式：

```matlab
M_physical = simplify(subs(M_a,pi_a,pi_param));
C_physical = simplify(subs(C_a,pi_a,pi_param));
G_physical = simplify(subs(G_a,pi_a,pi_param));
F_physical = simplify(subs(F_a,pi_a,pi_param));
```

只对新 `M_a` 代入质量和质心符号，不会自动替换其中的 `pi_a`。这一接口区别应在调用旧建模代码时明确处理。
