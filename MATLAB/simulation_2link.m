function result = simulation_2link(make_plots, settings)
%SIMULATION_2LINK 二连杆 AC 控制的 200 Hz 采样仿真。
%   simulation_2link                 % 运行并绘图
%   result = simulation_2link;        % 同时返回全部仿真数据
%   result = simulation_2link(false); % 不绘图，便于批量验证
%   result = simulation_2link(true,struct('t_end',20)); % 覆盖默认设置
%   settings 可设置 t_end、amplitude、omega、phase、Ks、Gamma。
%
%   AC_CON 是控制器；sim_model 是接收力矩、返回 q 和 dq 的仿真对象。
%   每 Ts=1/200=0.005 s 更新一次控制器，力矩在该区间内零阶保持，
%   对象由 sim_model 内部的 ode45 自适应龙格–库塔方法积分。
%   参数估计按采样自适应律更新：pi_hat(k+1)=pi_hat(k)+Ts*dpi_hat(k)。
%   这是 200 Hz 的采样控制，不是将控制器连续嵌入 ode45 的仿真。
%
%   默认 qd=[sin(t);sin(1.7*t+0.5)]，使用不同频率和相位改善激励。
%   dqd、ddqd 均为正弦轨迹的解析导数；默认仿真 120 秒。
%   若要复现原同频输入，可传 struct('omega',[1;1],'phase',[0;0])。
%   result 中每行对应一个采样时刻：
%     t                 (K+1)×1，时间，s
%     q、qd、e          (K+1)×2，位置、期望位置、位置误差，rad
%     dq、dqd、de       (K+1)×2，速度、期望速度、速度误差，rad/s
%     ddqd、s、tau      (K+1)×2，期望加速度、组合误差、控制力矩
%     pi_hat、dpi_hat   (K+1)×7，参数估计与参数估计导数
%     pi_true、pi_error、relative_pi_error 为仿真真值及估计误差诊断。
%   跟踪误差减小不代表参数已辨识；有限窗口的回归矩阵指标不是 PE 证明。
%   tau 最后一行仅用于记录终点控制器输出，未施加给对象。

if nargin < 1
    make_plots = true;
end
if nargin < 2 || isempty(settings)
    settings = struct();
end
if ~isstruct(settings) || ~isscalar(settings)
    error('simulation_2link:InvalidSettings', 'settings 必须为标量结构体。');
end
if ~((islogical(make_plots) || isnumeric(make_plots)) && ...
        isreal(make_plots) && isscalar(make_plots) && isfinite(make_plots) && ...
        (make_plots == 0 || make_plots == 1))
    error('simulation_2link:InvalidPlotFlag', 'make_plots 必须为 true 或 false。');
end
addpath(fileparts(mfilename('fullpath')));

%% 1. 仿真及正弦轨迹设置（在此修改）
fs = 200;                          % 控制采样频率，Hz
Ts = 1/fs;                         % 控制采样周期，s
t_end = setting_value(settings,'t_end',120); % 仿真时长，s
if ~isnumeric(t_end) || ~isreal(t_end) || ~isscalar(t_end) || ...
        ~isfinite(t_end) || t_end < Ts
    error('simulation_2link:InvalidDuration', 't_end 必须为不小于 Ts 的有限实数。');
end
step_count = round(t_end/Ts);
t = (0:step_count).'*Ts;
amplitude = reference_vector(setting_value(settings,'amplitude',[1;1]),'amplitude');
omega = reference_vector(setting_value(settings,'omega',[1;1.7]),'omega'); % rad/s
phase = reference_vector(setting_value(settings,'phase',[0;0.5]),'phase'); % rad
if any(omega <= 0)
    error('simulation_2link:InvalidFrequency', 'omega 必须为正角频率（rad/s）。');
end
angle = t*omega.' + phase.';
qd = sin(angle).*amplitude.';
dqd = cos(angle).*(amplitude.*omega).';
ddqd = -sin(angle).*(amplitude.*omega.^2).';

%% 2. 对象参数：与当前 sim_model 默认二连杆相同，显式传入避免旧模型残留
l1 = 0.45;
m = [23.902;3.88];
lc = [0.091;0.048];
Izz = [1.266;0.093];
fv = [2.288;0.175];
g = 9.81;
robot.MDH = [0,0,0,-pi/2;0,l1,0,0];
robot.m = m;
robot.c = [lc.';zeros(2,2)];
robot.Ic = zeros(3,3,2);
robot.Ic(:,:,1) = diag([0,0,Izz(1)]);
robot.Ic(:,:,2) = diag([0,0,Izz(2)]);
robot.fv = fv;
robot.gravity = [0;-g;0];

%% 3. 控制器及初值：Lambda 必须与当前 Y_a 的 3I 一致
ctrl.Lambda = 3*eye(2);
ctrl.Ks = setting_value(settings,'Ks',diag([20,10]));
ctrl.Gamma = setting_value(settings,'Gamma',diag([1,5,1,1,1,0.5,1]));
ctrl.g = g;
ctrl.adapt = true;

% 将表格标称参数按 AC_CON 的七列顺序组合；初值取 80% 制造估计偏差。
% 这是演示初始化。实际应用中用自己的标称估计，不在线读取对象参数。
pi_nominal = [fv(1); ...
    Izz(1)+m(1)*lc(1)^2+m(2)*l1^2; ...
    m(2)*l1*lc(2); ...
    m(1)*lc(1)+m(2)*l1; ...
    Izz(2)+m(2)*lc(2)^2; ...
    m(2)*lc(2);fv(2)];
pi_current = 0.8*pi_nominal;
pi_true = pi_nominal;               % 当前演示对象的真值，仅用于结果评估
q0 = zeros(2,1);
dq0 = zeros(2,1);
[dq_current,q_current] = sim_model(zeros(2,1),true,robot,q0,dq0);

%% 4. 200 Hz 控制循环：先记录/计算当前控制，再推进对象到下一采样时刻
sample_count = numel(t);
q = zeros(sample_count,2);
dq = zeros(sample_count,2);
e = zeros(sample_count,2);
de = zeros(sample_count,2);
s = zeros(sample_count,2);
tau = zeros(sample_count,2);
pi_hat = zeros(sample_count,7);
dpi_hat = zeros(sample_count,7);
excitation_start = max(0,t(end)-10);
gram = zeros(7);
fprintf('二连杆 AC 仿真：%.0f Hz，Ts=%.4f s，时长 %.1f s。\n',fs,Ts,t(end));
for k = 1:sample_count
    [tau_cmd,dpi_current,obs] = AC_CON(q_current,dq_current, ...
        qd(k,:).',dqd(k,:).',ddqd(k,:).',pi_current,ctrl);
    q(k,:) = q_current.';
    dq(k,:) = dq_current.';
    e(k,:) = obs.e.';
    de(k,:) = obs.de.';
    s(k,:) = obs.s.';
    tau(k,:) = tau_cmd.';
    pi_hat(k,:) = pi_current.';
    dpi_hat(k,:) = dpi_current.';
    if t(k) >= excitation_start && k < sample_count
        gram = gram+Ts*(obs.Y_a.'*obs.Y_a);
    end

    if k < sample_count
        % 本周期仅施加当前力矩；下一周期使用更新后的参数和对象反馈。
        pi_next = pi_current+Ts*dpi_current;
        if any(~isfinite(pi_next))
            error('simulation_2link:NonfiniteEstimate', ...
                't=%.6g s 时参数估计产生非有限值。',t(k));
        end
        [dq_current,q_current] = sim_model( ...
            tau_cmd,false,[],[],[],[t(k),t(k+1)]);
        pi_current = pi_next;
    end
end

%% 5. 返回数据与跟踪指标
rmse_q = sqrt(mean(e.^2,1));
rmse_dq = sqrt(mean(de.^2,1));
tail_index = t >= max(0,t(end)-5);
rmse_q_last5s = sqrt(mean(e(tail_index,:).^2,1));
pi_error = pi_hat-pi_true.';
relative_pi_error = pi_error./abs(pi_true).';
column_norm = sqrt(diag(gram));
if all(column_norm > 0)
    normalized_gram = gram./(column_norm*column_norm.');
    lambda_min = min(eig((normalized_gram+normalized_gram.')/2));
    condition_number = cond(normalized_gram);
else
    normalized_gram = zeros(7);
    lambda_min = 0;
    condition_number = Inf;
end
excitation = struct('window',[excitation_start,t(end)],'gram',gram, ...
    'normalized_gram',normalized_gram,'lambda_min',lambda_min, ...
    'condition_number',condition_number);
result = struct('t',t,'fs',fs,'Ts',Ts,'q',q,'dq',dq, ...
    'qd',qd,'dqd',dqd,'ddqd',ddqd,'e',e,'de',de,'s',s,'tau',tau, ...
    'pi_hat',pi_hat,'dpi_hat',dpi_hat,'pi_nominal',pi_nominal,'pi_true',pi_true, ...
    'pi_error',pi_error,'relative_pi_error',relative_pi_error, ...
    'excitation',excitation,'ctrl',ctrl, ...
    'reference',struct('amplitude',amplitude,'omega',omega,'phase',phase), ...
    'rmse_q',rmse_q,'rmse_dq',rmse_dq,'rmse_q_last5s',rmse_q_last5s);
fprintf('位置 RMSE：关节1 %.6f rad，关节2 %.6f rad。\n',rmse_q);
fprintf('最后5秒位置 RMSE：关节1 %.6f rad，关节2 %.6f rad。\n',rmse_q_last5s);
fprintf('速度 RMSE：关节1 %.6f rad/s，关节2 %.6f rad/s。\n',rmse_dq);
fprintf('参数编号       仿真真值        终点估计      相对误差(%%)\n');
for parameter = 1:7
    fprintf('  %d         %11.6f    %11.6f      %+9.3f\n', ...
        parameter,pi_true(parameter),pi_hat(end,parameter), ...
        100*relative_pi_error(end,parameter));
end
fprintf('最后10秒归一化回归 Gram：最小特征值 %.3e，条件数 %.3e。\n', ...
    lambda_min,condition_number);
fprintf('以上为有限时段诊断，不能据此保证参数严格收敛到真值。\n');
if make_plots
    plot_simulation(result);
end
end


function plot_simulation(result)
t = result.t;
figure('Name','二连杆 AC：位置、速度与误差','Color','w', ...
    'Position',[100,100,1100,780]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for joint = 1:2
    nexttile(joint);
    plot(t,result.qd(:,joint),'k--',t,result.q(:,joint),'b-','LineWidth',1.2);
    grid on; ylabel('q (rad)'); title(sprintf('关节 %d：位置跟踪',joint));
    legend('期望','实际','Location','best');
    nexttile(joint+2);
    plot(t,result.dqd(:,joint),'k--',t,result.dq(:,joint),'b-','LineWidth',1.2);
    grid on; ylabel('dq (rad/s)'); title(sprintf('关节 %d：速度跟踪',joint));
    legend('期望','实际','Location','best');
    nexttile(joint+4);
    plot(t,result.e(:,joint),'r-','LineWidth',1.2);
    grid on; xlabel('t (s)'); ylabel('e (rad)');
    title(sprintf('关节 %d：位置误差',joint));
end

figure('Name','二连杆 AC：力矩与组合误差','Color','w', ...
    'Position',[100,100,1000,660]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for joint = 1:2
    nexttile(joint);
    % stairs 对应 200 Hz 控制力矩的零阶保持。
    stairs(t,result.tau(:,joint),'b-','LineWidth',1.1);
    grid on; xlabel('t (s)'); ylabel('\tau (N m)');
    title(sprintf('关节 %d：控制力矩',joint));
    nexttile(joint+2);
    plot(t,result.s(:,joint),'r-','LineWidth',1.2);
    grid on; xlabel('t (s)'); ylabel('s (rad/s)');
    title(sprintf('关节 %d：s = de + 3e',joint));
end

figure('Name','二连杆 AC：七维参数估计','Color','w', ...
    'Position',[100,100,1100,850]);
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
parameter_names = {'f_{v1}','A','B','H_1','D','H_2','f_{v2}'};
for parameter = 1:7
    nexttile(parameter);
    plot(t,result.pi_hat(:,parameter),'b-','LineWidth',1.2);
    yline(result.pi_true(parameter),'k--','LineWidth',1.0);
    grid on; xlabel('t (s)'); ylabel('$\hat{\pi}$','Interpreter','latex');
    title(sprintf('参数 %d：%s',parameter,parameter_names{parameter}));
    legend('估计','仿真真值','Location','best');
end
nexttile(8);
plot(t,vecnorm(result.relative_pi_error,2,2),'r-','LineWidth',1.2);
grid on; xlabel('t (s)'); ylabel('相对参数误差范数');
title('七维参数误差（按真值归一化）');
end


function value = setting_value(settings, name, fallback)
if isfield(settings,name)
    value = settings.(name);
else
    value = fallback;
end
end


function value = reference_vector(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        numel(value) ~= 2 || any(~isfinite(value))
    error('simulation_2link:InvalidReference', '%s 必须为 2 维有限实数向量。',name);
end
value = double(value(:));
end
