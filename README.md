# AC普通自适应控制器MATLAB仿真
在MATLAB建立二连杆模型，进行自适应控制器仿真
## 前置代码
该代码的前置代码为：https://github.com/yzjjx/min_parameter  
前置代码实现了线性化模型：

$$
\tau = Y(q,\dot{q},\ddot{q})\pi
$$

## 后续实现
后续该代码主要实现机器人动力学的参数线性化模型，即

$$
\tau_a = Y_a\pi
$$

具体构建方式参考PPT文件：相关说明文档\普通AC构建.pptx

## 1 二连杆仿真模型建立
仿真模型建立：  
输入tau，输出q和dq  
具体公式： 

$$
\tau = M(q)\ddot{q} + C(q,\dot{q})\dot{q} + G(q) + F(\dot{q})
$$

主要使用RNEA计算，同时使用龙格库塔数值积分（ode45）方法  

主要步骤：  
1.输入当前 $q,\dot{q}$ 和控制力矩 $\tau$   
2.通过RNEA计算 $M(q)、C(q,\dot{q})\dot{q}、G(q) $   
3.求解加速度： $\ddot{q} = M^{-1}\left( \tau - C\dot{q} - G - F \right)$    
4.ode45在给定时间区间内积分运动微分方程  
5.输出并且保存下一时刻的 $q,\dot{q}$  

**注意**：龙格库塔（ode45）方法不是直接使用 $q_{next}=q+\dot{q} \Delta t$ 这种简单公式，而是通过多次计算微分方程，提高积分精度

具体文件参考：`MATLAB\sim_model.m`，该文件为m函数，用来配合控制器使用，可以在m函数内部设置N连杆参数模型  

主要龙格库塔代码：

```m
[t_solution,x_solution] = ode45( ...
    @(t,x) state_derivative(t,x,tau_cmd,candidate_robot), ...
    tspan,x_initial,options);
```

基本格式为：

```m
[t,x] = ode45(微分方程函数, 时间区间, 初始状态, 配置参数);
```

## 2 Y_a矩阵构建
核心方法：  
从已经提取出来的M、C和G、F矩阵构造tau_a的公式，然后对独立参数pi求雅可比得到Y_a  

前置代码Mat_out.m已经完成M、C和G矩阵的提取，后续通过摩擦力矩相减得到F，即代码：

```m
F = simplify(tau-tau_before_friction);
```

然后使用原始物理参数构造补偿力矩的公式，如下：

```m
e = mat_q-mat_qd;
de = mat_dq-mat_dqd;

v_r = mat_dqd-sym(Lambda)*e;
a_r = mat_ddqd-sym(Lambda)*de;

tau_a_physical = simplify(M*a_r+C*v_r+G+F);
```

然后将物理参数组合换成独立参数pi_a，在这里pi_a=pi，因此前置代码需要运行`MATLAB\generate_Y_pi.m`，同时，如果需要假如摩擦参数，需要先运行`MATLAB\add_friction.m`，再运行`MATLAB\generate_Y_pi.m`

最终公式即：
```m
tau_a = simplify(M_a*a_r + C_a*v_r + G_a + F_a);
Y_a = simplify(jacobian(tau_a, pi_a));
```

具体文件参考：`MATLAB\generate_Y_a.m`

## 3 控制器代码
得到Y_a的矩阵可以直接写控制器函数，可以参考PPT第四页，同时控制律如下公式：

$$
\begin{aligned}
\tau &= Y_a \hat{\pi} - K_s S \\
\dot{\hat{\pi}} &= -\Gamma Y_a^T S
\end{aligned}
$$

相关代码部分：

```m
tau_model = Y_a*pi_hat;
tau_feedback = -Ks*s;
tau_cmd = tau_model + tau_feedback;
```

具体文件参考：MATLAB\AC_CON.m  

注意：  
控制器需要时刻获取pi_hat的值，以观察参数是否完全收敛

## 4 最终结果
代码文件参考：MATLAB\simulation_2link.m  
### 20s仿真
仿真假设控制频率和发送频率均为200Hz，首先代码运行20s，查看参数收敛情况，运行结果如下图所示：

<div align="center">
    <img src="fig\20s_sim\pi_hat_20s.png">
    <br>
    图1：参数收敛情况
</div>


<div align="center">
    <img src="fig\20s_sim\pos_err_20s.png">
    <br>
    图2：位置误差
</div>

### 120s仿真

<div align="center">
    <img src="fig\120s_sim\pi_hat_120s.png">
    <br>
    图3：参数收敛情况
</div>


<div align="center">
    <img src="fig\120s_sim\pos_err_120s.png">
    <br>
    图4：位置误差
</div>
