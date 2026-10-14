仿真模型建立：  
输入tau，输出q和dq  
具体公式：  
$$

$$

然后对加速度两次积分，即可得到速度和位置  
具体文件参考：MATLAB\sim_model.m，该文件为m函数，用来配合控制器使用，可以在m函数内部设置N连杆参数模型  
使用牛顿欧拉计算机器人动力学，再求解关节加速度，最后通过ode45积分得到机器人下一时刻的位置和速度
主要代码：  

## Y_a矩阵构建

```m
tau_a = simplify(M_a*a_r + C_a*v_r + G_a + F_a);
Y_a = simplify(jacobian(tau_a, pi_a));
```