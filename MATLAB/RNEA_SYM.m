% 设机器人标准动力学方程为：tau = M*ddq+C*dq+G

% 定义符号变量
syms g real

% 定义关节数量-------------------------------------------------------------------->调整模型需要修改
N = 2;

% 关节旋转角矩阵，使用sym函数，根据关节数量N自动生成符号变量列向量
% 以二关节为例，mat_q=[q_1;q_2],下面定义关节位置、关节速度、关节加速度
mat_q = sym('q_',[N,1]);
mat_dq = sym('dq_', [N, 1]);
mat_ddq = sym('ddq_', [N, 1]);

% 连杆长度矩阵
mat_L = sym('L_',[N,1]);

% 输入：MDH表-------------------------------------------------------------------->调整模型需要修改
% Mat_MDH = [
%     0 0         0   mat_q(1);
%     0 mat_L(1)  0   mat_q(2)];

% mat_q 现在代表论文定义的关节角（初始0度在上）
Mat_MDH = [
    0 0         0   mat_q(1) - sym(pi)/2;
    0 mat_L(1)  0   mat_q(2)
];

% MDH建模，设轴均为局部+z方向,在建模时一定注意建模方向===============================>建模过程需要注意
% 所有关节均为旋转关节，旋转轴沿各自局部 +Z 轴
% mat_S 为 6×N，每一列对应一个关节的 S_i
mat_S = sym(zeros(6, N));
mat_S(3, :) = 1;

% 参考书本公式5.16，计算相对空间速度v_j
% v_{Ji}=S_i\dot{q}_i
% mat_vJ 为 6×N，每一列对应一个关节的 v_Ji
mat_vJ = mat_S * diag(mat_dq);

% 空间变换
% 齐次变换矩阵函数，根据MDH得到齐次变换矩阵T，T^1代表 ^1T_2 这个齐次变换矩阵
T = make_T(Mat_MDH);% 自定义函数使用
% 元胞数组存储矩阵
mat_T = cell(N,1);
for i = 1:N
    mat_T{i} = T(:,:,i);
end
% 调用示例：
% mat_T{1} 对应 T(:,:,1)
% mat_T{2} 对应 T(:,:,2)

% 生成6*6空间变换矩阵（参考书本：刚体动力学）
% 注意：这里的X_1对应公式的^1X_{\lambda (1)}
X = make_X(T);% 自定义函数使用
mat_X = cell(N,1);
for i = 1:N
    mat_X{i} = X(:,:,i);
end

% 初始化机器人基坐标系
% 机器人基座固定，v_0=0
v_0 = [0;0;0;0;0;0];

% 设定世界重力加速度，目前为沿y方向向下----------------------------------------------->调整模型需要修改
mat_g = [0;-g;0];

% 计算机器人基座虚拟加速度
a_0 = -[0;0;0;mat_g];

% ======================================================================
% 前向递推：从基座到末端
% 运动学递推（速度与加速度）
mat_v = cell(N,1); % 保存各关节的空间速度，mat_v{i} 对应 v_i
mat_a = cell(N,1); % 保存各关节的空间加速度，mat_a{i} 对应 a_i
for i = 1:N
    % ============计算速度================
    % 参考公式5.18：v_i=^iX_{\lambda (i)}v_{\lambda (i)}+v_{Ji}
    % 因为单独设置基座的速度和加速度，所以这里需要判断，不然直接循环就行
    if i == 1
        v_i = mat_X{i}*v_0 + mat_vJ(:,i);
    else
        v_i = mat_X{i}*mat_v{i-1} + mat_vJ(:,i);
    end
    mat_v{i} = simplify(v_i);
    % 自动生成工作区变量 v_1、v_2、...、v_N，供直接查看
    eval(sprintf('v_%d = mat_v{%d};', i, i));

    % ============计算加速度===============
    % 参考公式5.19：a_i = ^iX_{\lambda (i)}a_{\lambda (i)}+...
    % S_i\ddot{q}_i+c_{Ji}+v_i \times  v_{Ji}(接上)
    % 当前关节的 S_i 为常量，因此 c_Ji = 0，不加入计算
    % 同时，因为是6*1矩阵叉乘6*1矩阵，该计算法则的定义在书本2.33公式
    % 这里直接使用函数
    if i == 1
        a_i = mat_X{i}*a_0 + mat_S(:,i)*mat_ddq(i) + ...
              cross_6x1(mat_v{i}, mat_vJ(:,i));
    else
        a_i = mat_X{i}*mat_a{i-1} + mat_S(:,i)*mat_ddq(i) + ...
              cross_6x1(mat_v{i}, mat_vJ(:,i));
    end
    mat_a{i} = simplify(a_i);
    % 自动生成工作区变量 a_1、a_2、...、a_N，供直接查看
    eval(sprintf('a_%d = mat_a{%d};', i, i));
end

% 动力学递推(力与力矩)
% 建立空间惯性矩阵
% 在建模过程中，坐标系原点设置在关节，与质心不重合，使用公式2.63
% 公式2.63 latex：
% I_O=\begin{bmatrix}
%  \bar{I}_C+m\mathbf{c}\times\mathbf{c}\times^T    & m\mathbf{c}\times \\
%  m\mathbf{c}\times^T  &m\mathbf{1}
% \end{bmatrix}

% 其中c=[cx;cy;cz]为质心相对于连杆坐标系的位置向量
% 每一列为一个连杆的质心位置 [cx_i; cy_i; cz_i]
mat_c = [
    sym('cx_', [1, N]);
    sym('cy_', [1, N]);
    sym('cz_', [1, N])
];

% \bar{I}_C表示刚体相对于质心C的转动惯量矩阵
I_xx = sym('I_xx', [1, N]);
I_xy = sym('I_xy', [1, N]);
I_xz = sym('I_xz', [1, N]);
I_yy = sym('I_yy', [1, N]);
I_yz = sym('I_yz', [1, N]);
I_zz = sym('I_zz', [1, N]);

mat_Ic = sym(zeros(3, 3, N));
for i = 1:N
    mat_Ic(:,:,i) = [
        I_xx(i), I_xy(i), I_xz(i);
        I_xy(i), I_yy(i), I_yz(i);
        I_xz(i), I_yz(i), I_zz(i)
    ];
end

% m\mathbf{1}表示质量m乘以三阶单位矩阵
% 连杆质量列向量
mat_m = sym('m_', [N, 1]);
% 每个连杆的质量乘以三阶单位矩阵
mat_mI = sym(zeros(3, 3, N));
for i = 1:N
    mat_mI(:,:,i) = mat_m(i) * eye(3);
end

% 完成叉乘矩阵函数c*，建立六维惯量矩阵
% mat_Ic 和 mat_c 均在对应连杆坐标系的坐标轴方向下表达
mat_I = sym(zeros(6, 6, N));
for i = 1:N
    c_cross = cross_3x1(mat_c(:,i));
    m_i = mat_m(i);

    mat_I(:,:,i) = simplify([ ...
        mat_Ic(:,:,i) + m_i*c_cross*c_cross.', m_i*c_cross; ...
        m_i*c_cross.', mat_mI(:,:,i)]);

    % 自动生成工作区变量 I_1、I_2、...、I_N，供直接查看
    eval(sprintf('I_%d = mat_I(:,:,%d);', i, i));
end

mat_fB = cell(N,1); % 保存各连杆的本体空间力，排列为 [力矩; 力]
for i = 1:N
    %===========================计算本体空间力f_i^B======================
    % f_i^B = I_i*a_i + v_i ×* (I_i*v_i)
    % 注：空间力矩阵为6*1矩阵，前三行绕XYZ轴的力矩，后三行为XYZ方向的力
    fB_i = mat_I(:,:,i)*mat_a{i} + ...
           cross_star(mat_v{i}, mat_I(:,:,i)*mat_v{i});
    mat_fB{i} = simplify(fB_i);
    % 自动生成工作区变量 fB_1、fB_2、...、fB_N，供直接查看
    eval(sprintf('fB_%d = mat_fB{%d};', i, i));
end

% ======================================================================
% 反向递推：从末端到基座
% 串联机器人：连杆 i 的子连杆为 i+1，当前外部空间力 f_i^x = 0
tau = sym(zeros(N,1));
mat_f = cell(N,1); % 保存含下游载荷的空间力，排列为 [力矩; 力]
for i = N:-1:1
    if i == N
        % 最后连杆无子连杆，f_i=fB_i
        mat_f{i} = mat_fB{i};
    else
        % 公式5.20：f_i = f_i^B + X_{i+1}.'*f_{i+1}
        % mat_X{i+1} 为父系 i 到子系 i+1 的空间运动变换
        mat_f{i} = simplify(mat_fB{i} + mat_X{i+1}.'*mat_f{i+1});
    end

    % 公式5.11：tau_i = S_i.'*f_i，取运动子空间的第 i 列
    tau(i) = simplify(mat_S(:,i).'*mat_f{i});
    % 自动生成 f_1、...、f_N 和 tau_1、...、tau_N，供直接查看
    eval(sprintf('f_%d = mat_f{%d};', i, i));
    eval(sprintf('tau_%d = tau(%d);', i, i));
end
