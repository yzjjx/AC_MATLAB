function [dq_out, q_out] = sim_model(tau_cmd, reset, robot_params, q_init, dq_init, tspan)
%SIM_MODEL N 个旋转关节串联机器人的正向动力学，采用龙格–库塔积分。
%   [dq_out,q_out] = sim_model(tau_cmd,false,[],[],[],[t0 tf])
%   输入 tau_cmd：N×1 关节力矩，N*m；输出 dq_out：rad/s，q_out：rad。
%   每次普通调用由 ode45 在 tspan=[t0 tf]（秒）内自适应积分。
%   不设置固定采样周期；输出为 tf 时刻的状态。
%   函数保存内部状态，不依赖 Symbolic Math Toolbox，不等待实际时钟。
%
%   初始化/重置（只返回初始状态，不积分、不施加本次 tau_cmd）：
%     [dq,q] = sim_model([0;0],true);  % 默认表格中的二连杆，零初值
%     [dq,q] = sim_model(zeros(N,1),true,robot,q0,dq0);
%   恒定力矩（在整个积分区间内保持不变）：
%     [dq,q] = sim_model(tau_cmd,false,[],[],[],[0 1]);
%   连续时间控制（每个积分阶段重新计算力矩）：
%     torque_fun = @(t,q,dq) controller(t,q,dq);
%     [dq,q] = sim_model(torque_fun,false,[],[],[],[0 1]);
%   分段仿真时，下一段从上一段的 tf 开始，使用保存的终态。
%   torque_fun 须返回 N×1 有限实数力矩；不得在回调中累计内部状态，
%   因为 ode45 可能重复计算或拒绝试探步。
%
%   robot_params 为标量结构体，换 N 或换模型时通过 reset=true 提交：
%     MDH       N×4，每行 [alpha,a,d,theta_offset]（固定部分）
%               theta_i = theta_offset_i + q_i，变换与 make_T.m 一致
%     m         N×1 连杆质量，kg
%     c         3×N 质心在各连杆坐标系中的位置，m
%     Ic        3×3×N 关于质心的惯量，在对应连杆坐标系表达，kg*m^2
%     fv        N×1 黏性摩擦系数，N*m*s/rad
%     gravity   3×1 世界重力加速度，m/s^2
%   未提供 robot_params 时，首次调用使用下方显式表格参数。
%   重置已有模型而未提供 robot_params 时，保留该模型并重置状态。
%   表格只有两根杆的数据；N>2 时须提供完整 N 杆参数，不能只改力矩长度。
%
%   动力学（沿用 RNEA_SYM 的固定基座、局部 +z 旋转轴约定）：
%     tau = M(q)*ddq + C(q,dq)*dq + G(q) + fv.*dq
%     ddq = M(q) \ (tau_cmd - C_dq - G - F)
%   一阶状态 x=[q;dq]，dx=[dq;ddq]；使用 ode45 的显式 RK(4,5)。
%   每个积分阶段都按当前 q、dq 重新计算 M、C_dq、G 和 F。
%   q_init、dq_init 只在首次初始化或 reset=true 时使用，默认均为零。
%   新实验可 reset=true；clear sim_model 可清除状态及保存的模型。

%% 1. 输入与保存的状态

% yzj_add_annotation：首先持久化局部变量，persistent声明的便里昂可以再下一次调用同一个函数的时候继续使用
% yzj_add_annotation：q_state和dq_state分别保存了上一次调用sim_model函数时的关节位置和速度状态，robot_state保存了机器人模型参数，包括空间惯量等
persistent q_state dq_state robot_state

% yzj_add_annotation：判断
if nargin < 2 || isempty(reset)
    reset = false;
end
if ~((islogical(reset) || isnumeric(reset)) && isreal(reset) && ...
        isscalar(reset) && isfinite(reset) && (reset == 0 || reset == 1))
    error('sim_model:InvalidReset', 'reset 必须是 true 或 false。');
end
reset = logical(reset);
if nargin < 3
    robot_params = [];
end
if nargin < 4
    q_init = [];
end
if nargin < 5
    dq_init = [];
end
if nargin < 6
    tspan = [];
end
initializing = reset || isempty(robot_state);
if ~initializing && (~isempty(robot_params) || ~isempty(q_init) || ~isempty(dq_init))
    error('sim_model:ResetRequired', ...
        '更换模型或初值时请设置 reset=true；普通调用的模型和初值传 []。');
end

%% 2. 准备 N 杆参数；空间惯量只在初始化/重置时计算
if initializing
    if ~isempty(robot_params)
        candidate_robot = prepare_robot(robot_params);
    elseif ~isempty(robot_state)
        candidate_robot = robot_state;
    else
        candidate_robot = prepare_robot(default_robot_parameters());
    end
    N = candidate_robot.N;
    if isempty(q_init)
        q_init = zeros(N,1);
    end
    if isempty(dq_init)
        dq_init = zeros(N,1);
    end
    q_current = state_vector(q_init, N, 'q_init');
    dq_current = state_vector(dq_init, N, 'dq_init');
else
    candidate_robot = robot_state;
    N = candidate_robot.N;
    q_current = q_state;
    dq_current = dq_state;
end
if ~isa(tau_cmd, 'function_handle')
    tau_cmd = state_vector(tau_cmd, N, 'tau_cmd');
end

if reset
    robot_state = candidate_robot;
    q_state = q_current;
    dq_state = dq_current;
    dq_out = dq_state;
    q_out = q_state;
    return;
end

%% 3. 自适应龙格–库塔积分；tspan 仅指定起止时间，不指定内部步长
if ~isnumeric(tspan) || ~isreal(tspan) || ~isvector(tspan) || ...
        numel(tspan) ~= 2 || any(~isfinite(tspan)) || tspan(2) <= tspan(1)
    error('sim_model:InvalidTimeSpan', ...
        '普通调用必须提供 tspan=[t0 tf]，且 t0、tf 有限、tf>t0。');
end
tspan = double(tspan(:).');
x_initial = [q_current;dq_current];
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t_solution,x_solution] = ode45( ...
    @(t,x) state_derivative(t,x,tau_cmd,candidate_robot), ...
    tspan,x_initial,options);
if t_solution(end) ~= tspan(2)
    error('sim_model:IntegrationIncomplete', ...
        '积分未到达 tf，模型状态未更新。');
end
x_next = x_solution(end,:).';
if any(~isfinite(x_next))
    error('sim_model:NonfiniteState', '计算产生非有限状态，模型状态未更新。');
end

% 全部计算成功后再提交状态；输出是 tf 时刻的反馈。
robot_state = candidate_robot;
q_state = x_next(1:N);
dq_state = x_next(N+1:2*N);
dq_out = dq_state;
q_out = q_state;
end


function dx = state_derivative(t, x, tau_cmd, robot)
N = robot.N;
if any(~isfinite(x))
    error('sim_model:NonfiniteState', '计算产生非有限状态，模型状态未更新。');
end
q = x(1:N);
dq = x(N+1:2*N);
if isa(tau_cmd, 'function_handle')
    tau = state_vector(tau_cmd(t,q,dq), N, 'torque_fun 返回值');
else
    tau = tau_cmd;
end

% 每个 RK 阶段使用该阶段的状态重建动力学。
X = motion_transforms(robot.MDH, q);

% 提取 tau 对 ddq 的系数。零速度、零重力下，N 列单位加速度
% 同时通过 RNEA 传播，相当于并行计算 tau(q,0,e_j,0) 的所有列。
M = acceleration_coefficients(X, robot.spatial_I, N);

% 无需显式存储 C 矩阵：分别用 RNEA 得到 C(q,dq)*dq 与 G(q)。
C_dq = rnea_bias(X, robot.spatial_I, dq, zeros(3,1), N);
G = rnea_bias(X, robot.spatial_I, zeros(N,1), robot.gravity, N);
F = robot.fv .* dq;

[~, mass_flag] = chol(M);
if mass_flag ~= 0 || any(~isfinite(M), 'all')
    error('sim_model:InvalidMassMatrix', ...
        '当前惯性矩阵不是有限正定矩阵，请检查 MDH、质心和惯量。');
end

ddq = M \ (tau - C_dq - G - F);
dx = [dq;ddq];
if any(~isfinite(dx))
    error('sim_model:NonfiniteState', '计算产生非有限状态，模型状态未更新。');
end

end


function robot = default_robot_parameters()
% 表格参数显式列出：默认二连杆。替换其他机器人时传入 robot_params。
l1 = 0.45;                          % 连杆 1 长度，m
l2 = 0.45;                          % 连杆 2 长度，m
lc1 = 0.091;                        % 连杆 1 质心距离，m
lc2 = 0.048;                        % 连杆 2 质心距离，m
m1 = 23.902;                        % 连杆 1 质量，kg
m2 = 3.88;                          % 连杆 2 质量，kg
I1 = 1.266;                         % 连杆 1 质心 z 轴惯量，kg*m^2
I2 = 0.093;                         % 连杆 2 质心 z 轴惯量，kg*m^2
g = 9.81;                          % 重力加速度，m/s^2
fv1 = 2.288;                        % 关节 1 黏性摩擦，N*m*s/rad
fv2 = 0.175;                        % 关节 2 黏性摩擦，N*m*s/rad

% 与当前 RNEA_SYM 的第一关节 -pi/2 偏移和负 y 重力一致。
robot.MDH = [0, 0,  0, -pi/2; ...
             0, l1, 0,  0];
robot.m = [m1;m2];
robot.c = [lc1,lc2;0,0;0,0];
robot.Ic = zeros(3,3,2);
robot.Ic(:,:,1) = diag([0,0,I1]);
robot.Ic(:,:,2) = diag([0,0,I2]);
robot.fv = [fv1;fv2];
robot.gravity = [0;-g;0];
robot.link_lengths = [l1;l2];
% 表格只给 Izz，默认按 data_verify 的平面模型设置其他惯量分量为零。
% 非平面 N 杆模型须传入完整 Ic。末杆长度 l2 已记录；在当前无末端
% 外力模型中，它不独立进入关节力矩，其质量分布由 lc2、I2 描述。
end


function robot = prepare_robot(robot)
if ~isstruct(robot) || ~isscalar(robot)
    error('sim_model:InvalidRobot', 'robot_params 必须是标量结构体。');
end
required = {'MDH','m','c','Ic','fv','gravity'};
for k = 1:numel(required)
    if ~isfield(robot, required{k})
        error('sim_model:MissingRobotParameter', 'robot_params 缺少 %s。', required{k});
    end
    value = robot.(required{k});
    if ~isnumeric(value) || ~isreal(value) || isempty(value) || ...
            any(~isfinite(value), 'all')
        error('sim_model:InvalidRobotParameter', '%s 必须是有限实数数组。', required{k});
    end
    robot.(required{k}) = double(value);
end
if ~isvector(robot.m) || any(robot.m <= 0)
    error('sim_model:InvalidMass', 'm 必须是严格正值质量向量。');
end
robot.m = robot.m(:);
N = numel(robot.m);
if ~isequal(size(robot.MDH), [N,4]) || ~isequal(size(robot.c), [3,N]) || ...
        size(robot.Ic,1) ~= 3 || size(robot.Ic,2) ~= 3 || ...
        size(robot.Ic,3) ~= N || ndims(robot.Ic) > 3
    error('sim_model:RobotDimensionMismatch', ...
        'm 为 N×1，MDH 须为 N×4，c 须为 3×N，Ic 须为 3×3×N。');
end
robot.fv = state_vector(robot.fv, N, 'fv');
if any(robot.fv < 0)
    error('sim_model:InvalidFriction', '黏性摩擦系数 fv 必须非负。');
end
robot.gravity = state_vector(robot.gravity, 3, 'gravity');
robot.N = N;
robot.spatial_I = zeros(6,6,N);
for i = 1:N
    Ic = robot.Ic(:,:,i);
    tolerance = 1e-10 * max(1,norm(Ic,'fro'));
    if norm(Ic-Ic.','fro') > tolerance || ...
            min(eig((Ic+Ic.')/2)) < -tolerance
        error('sim_model:InvalidInertia', 'Ic(:,:,%d) 必须对称且半正定。', i);
    end
    Ic = (Ic+Ic.')/2;
    cx = skew3(robot.c(:,i));
    mass = robot.m(i);
    robot.spatial_I(:,:,i) = [Ic+mass*cx*cx.', mass*cx; ...
                             mass*cx.', mass*eye(3)];
end
end


function v = state_vector(v, N, name)
if ~isnumeric(v) || ~isreal(v) || ~isvector(v) || numel(v) ~= N || ...
        any(~isfinite(v))
    error('sim_model:InvalidVector', '%s 必须是长度为 %d 的有限实数向量。', name, N);
end
v = double(v(:));
end


function X = motion_transforms(MDH, q)
N = numel(q);
X = zeros(6,6,N);
for i = 1:N
    alpha = MDH(i,1);
    theta = MDH(i,4) + q(i);
    ca = cos(alpha); sa = sin(alpha);
    ct = cos(theta); st = sin(theta);
    R = [ct,-st,0;ca*st,ca*ct,-sa;sa*st,sa*ct,ca];
    p = [MDH(i,2);-MDH(i,3)*sa;MDH(i,3)*ca];
    Rt = R.';
    X(:,:,i) = [Rt,zeros(3);-Rt*skew3(p),Rt];
end
end

% yzj_add_annotation：计算质量矩阵的函数
function M = acceleration_coefficients(X, spatial_I, N)
% RNEA 在零速度、零重力下对 N 组单位加速度同时计算。
% acceleration(:,j) 是 ddq=e_j 时当前连杆的空间加速度。
acceleration = zeros(6,N);
forces = zeros(6,N,N);
for i = 1:N
    acceleration = X(:,:,i)*acceleration;
    acceleration(3,i) = acceleration(3,i) + 1;
    forces(:,:,i) = spatial_I(:,:,i)*acceleration;
end
M = zeros(N);
for i = N:-1:1
    M(i,:) = forces(3,:,i);          % S_i.'*f_i，S_i 的第 3 分量为 1
    if i > 1
        forces(:,:,i-1) = forces(:,:,i-1) + X(:,:,i).'*forces(:,:,i);
    end
end
M = (M+M.')/2;                      % 消除浮点舍入带来的微小非对称
end


function tau_bias = rnea_bias(X, spatial_I, dq, gravity, N)
% RNEA_SYM 的数值版本，此处 ddq=0。
% gravity=0 时得到 C(q,dq)*dq；dq=0 时得到 G(q)。
velocity = zeros(6,1);
acceleration = [zeros(3,1);-gravity];
forces = zeros(6,N);
for i = 1:N
    vJ = [0;0;dq(i);0;0;0];
    velocity = X(:,:,i)*velocity + vJ;
    omega_cross = skew3(velocity(1:3));
    linear_cross = skew3(velocity(4:6));
    motion_cross = [omega_cross,zeros(3);linear_cross,omega_cross];
    acceleration = X(:,:,i)*acceleration + motion_cross*vJ;
    force_cross = [omega_cross,linear_cross;zeros(3),omega_cross];
    forces(:,i) = spatial_I(:,:,i)*acceleration ...
                  + force_cross*(spatial_I(:,:,i)*velocity);
end
tau_bias = zeros(N,1);
for i = N:-1:1
    tau_bias(i) = forces(3,i);
    if i > 1
        forces(:,i-1) = forces(:,i-1) + X(:,:,i).'*forces(:,i);
    end
end
end


function S = skew3(v)
S = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
end
