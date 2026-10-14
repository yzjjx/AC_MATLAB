function [tau_cmd, dpi_hat, obs] = AC_CON(q, dq, qd, dqd, ddqd, pi_hat, ctrl)
%AC_CON 二连杆七参数普通自适应控制器（Slotine–Li 型）。
%   [tau_cmd,dpi_hat,obs] = AC_CON(q,dq,qd,dqd,ddqd,pi_hat)
%   [tau_cmd,dpi_hat,obs] = AC_CON(q,dq,qd,dqd,ddqd,pi_hat,ctrl)
%
%   输入（行/列向量均可）：
%     q、dq          2 维实际位置、速度，rad、rad/s
%     qd、dqd、ddqd  2 维期望位置、速度、加速度，rad、rad/s、rad/s^2
%     pi_hat         7 维当前参数估计，必须与下方 Y_a 的列顺序一致
%     ctrl           可选标量结构体，未填写的字段使用以下默认值：
%       Ks           diag([8,4])，误差反馈增益，2×2 对称正定矩阵
%       Gamma        diag([0.01,0.3,0.02,0.02,0.2,0.01,0.01])，7×7 正定矩阵
%       g            9.81，已知重力加速度大小，m/s^2
%       adapt        true；false 时冻结参数估计（dpi_hat=0）
%       Lambda       3*eye(2)，已固化在当前 Y_a 中，不可单独修改
%       Ks、Gamma、Lambda 也可传正标量，表示该标量乘单位矩阵。
%
%   输出：tau_cmd 为 2×1 关节力矩（N*m），dpi_hat 为 7×1 参数估计导数，
%   obs 包含 e、de、s、v_r、a_r、Y_a、tau_model 和 tau_feedback。
%
%   本函数数值实现 generate_Y_a.m 输出的当前 2×7 矩阵，无符号运算。
%   参数顺序（I1、I2 均为质心处的 z 轴惯量）：
%     [fv1;
%      I1 + m1*lc1^2 + m2*l1^2;
%      m2*l1*lc2;
%      m1*lc1 + m2*l1;
%      I2 + m2*lc2^2;
%      m2*lc2;
%      fv2]
%   修改前置参数基或 Lambda 后，须同步重新生成本函数的 Y_a。
%
%   控制律及自适应律：
%     e = q-qd; de = dq-dqd; s = de+Lambda*e;
%     tau_cmd = Y_a*pi_hat - Ks*s;
%     dpi_hat = -Gamma*Y_a.'*s;
%
%   调用示例（pi_hat0 应由标称参数按上述顺序构造，也可从零估计开始）：
%     pi_hat0 = zeros(7,1);
%     [tau_cmd,dpi_hat,obs] = AC_CON( ...
%         [0.1;-0.1],[0;0],[0;0],[0;0],[0;0],pi_hat0);
%
%   连续仿真时将 pi_hat 纳入 ode45 状态：x=[q;dq;pi_hat]（11×1）。
%   在 ODE 右端调用本函数，使用对象动力学求 ddq，返回 dx=[dq;ddq;dpi_hat]。
%   dpi_hat 是变化率，不是下一时刻估计；本函数不保存状态或自行积分，
%   可被龙格–库塔的多个阶段反复调用。勿在右端调用保存状态的 sim_model。

%% 1. 输入与控制增益
q = finite_vector(q, 2, 'q');
dq = finite_vector(dq, 2, 'dq');
qd = finite_vector(qd, 2, 'qd');
dqd = finite_vector(dqd, 2, 'dqd');
ddqd = finite_vector(ddqd, 2, 'ddqd');
pi_hat = finite_vector(pi_hat, 7, 'pi_hat');

if nargin < 7 || isempty(ctrl)
    ctrl = struct();
end
if ~isstruct(ctrl) || ~isscalar(ctrl)
    error('AC_CON:InvalidConfig', 'ctrl 必须为标量结构体。');
end

Ks = positive_gain(option_value(ctrl, 'Ks', diag([8,4])), 2, 'Ks');
Gamma = positive_gain(option_value(ctrl, 'Gamma', ...
    diag([0.01,0.3,0.02,0.02,0.2,0.01,0.01])), 7, 'Gamma');
Lambda = positive_gain(option_value(ctrl, 'Lambda', 3*eye(2)), 2, 'Lambda');

if ~isequal(Lambda, 3*eye(2))
    error('AC_CON:LambdaMismatch', ...
        '当前 Y_a 已代入 Lambda=3*eye(2)，修改 Lambda 须同步重新生成 Y_a。');
end

g = option_value(ctrl, 'g', 9.81);

if ~isnumeric(g) || ~isreal(g) || ~isscalar(g) || ~isfinite(g) || g < 0
    error('AC_CON:InvalidGravity', 'g 必须为有限非负实数标量。');
end

g = double(g);
adapt = option_value(ctrl, 'adapt', true);

if ~((islogical(adapt) || isnumeric(adapt)) && isreal(adapt) && ...
        isscalar(adapt) && isfinite(adapt) && (adapt == 0 || adapt == 1))
    error('AC_CON:InvalidAdapt', 'adapt 必须为 true 或 false。');
end

%% 2. 跟踪误差与参考速度/加速度；与前置 generate_Y_a.m 一致
e = q - qd;
de = dq - dqd;
v_r = dqd - Lambda*e;
a_r = ddqd - Lambda*de;
s = de + Lambda*e;

%% 3. 数值 Y_a（严格保留用户给定的七列顺序）
% 使用 v_r、a_r 简写已展开的误差表达式，不再重复代入 Lambda。
c2 = cos(q(2));
s2 = sin(q(2));
sum_ar = a_r(1) + a_r(2);
coupling_1 = 2*c2*a_r(1) + c2*a_r(2) ...
    - dq(2)*s2*v_r(1) - (dq(1)+dq(2))*s2*v_r(2);
coupling_2 = c2*a_r(1) + dq(1)*s2*v_r(1);
gravity_1 = g*sin(q(1));
gravity_12 = g*sin(q(1)+q(2));

Y_a = [dq(1), a_r(1), coupling_1, gravity_1, sum_ar, gravity_12, 0; ...
       0,     0,      coupling_2, 0,         sum_ar, gravity_12, dq(2)];

%% 4. 力矩与参数估计导数；由外部 ODE 求解器统一积分
tau_model = Y_a*pi_hat;
tau_feedback = -Ks*s;
tau_cmd = tau_model + tau_feedback;
if adapt
    dpi_hat = -Gamma*Y_a.'*s;
else
    dpi_hat = zeros(7,1);
end
if any(~isfinite([e;de;v_r;a_r;s;Y_a(:);tau_cmd;dpi_hat]))
    error('AC_CON:NonfiniteOutput', '控制计算产生非有限结果，请检查输入与增益。');
end

if nargout > 2
    obs = struct('e',e,'de',de,'s',s,'v_r',v_r,'a_r',a_r, ...
        'Y_a',Y_a,'tau_model',tau_model,'tau_feedback',tau_feedback);
end
end


function v = finite_vector(v, count, name)
if ~isnumeric(v) || ~isreal(v) || ~isvector(v) || numel(v) ~= count || ...
        any(~isfinite(v))
    error('AC_CON:InvalidVector', ...
        '%s 必须为长度 %d 的有限实数向量。', name, count);
end
v = double(v(:));
end


function value = option_value(ctrl, name, fallback)
if isfield(ctrl, name)
    value = ctrl.(name);
else
    value = fallback;
end
end


function gain = positive_gain(gain, count, name)
if ~isnumeric(gain) || ~isreal(gain) || isempty(gain) || ...
        any(~isfinite(gain), 'all')
    error('AC_CON:InvalidGain', '%s 必须为有限实数正定增益。', name);
end
gain = double(gain);
if isscalar(gain)
    gain = gain*eye(count);
end
if ~isequal(size(gain), [count,count]) || ...
        norm(gain-gain.','fro') > 1e-12*max(1,norm(gain,'fro'))
    error('AC_CON:InvalidGain', '%s 必须为 %d×%d 对称正定矩阵。', ...
        name, count, count);
end
gain = (gain+gain.')/2;
[~,flag] = chol(gain);
if flag ~= 0
    error('AC_CON:InvalidGain', '%s 必须为对称正定矩阵。', name);
end
end
