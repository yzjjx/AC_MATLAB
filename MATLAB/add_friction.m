% 读取工作区已有 tau，加入粘性摩擦后仍输出 tau。
% 运行顺序：RNEA_SYM; data_verify; add_friction;
% 不清空工作区，也不重新计算刚体动力学。
% 修改参数或再次运行本脚本前，应重新运行前置脚本，避免重复叠加摩擦。
% 含摩擦的驱动力矩：tau = M*ddq + C*dq + G + Fv.*dq。
% 摩擦对关节的作用力矩为 -Fv.*dq，因此所需驱动力矩中取正号。

required_vars = {'tau', 'N'};
for k = 1:numel(required_vars)
    if ~exist(required_vars{k}, 'var')
        error('add_friction:MissingWorkspaceVariable', ...
            '工作区缺少变量 %s，请先运行 RNEA_SYM.m 和 data_verify.m。', required_vars{k});
    end
end
if ~(isa(tau, 'sym') || isnumeric(tau)) || ~isequal(size(tau), [N, 1])
    error('add_friction:InvalidTau', '工作区 tau 必须是 N×1 符号或数值列向量。');
end

%% 关节速度：与前置脚本对 tau 使用的速度保持一致
if exist('mat_dq_value', 'var')
    dq_friction = mat_dq_value; % data_verify 中实际代入的符号或数值速度
elseif exist('mat_dq', 'var')
    dq_friction = mat_dq;       % 仅运行 RNEA_SYM 时使用原始符号速度
else
    error('add_friction:MissingVelocity', '工作区需要 mat_dq_value 或 mat_dq。');
end
if ~(isa(dq_friction, 'sym') || isnumeric(dq_friction)) || ...
        ~isequal(size(dq_friction), [N, 1])
    error('add_friction:InvalidVelocity', '关节速度必须是 N×1 符号或数值列向量。');
end

%% 粘性摩擦：每个关节一个系数，旋转关节单位为 N*m*s/rad
mat_Fv_value = sym('Fv_', [N, 1], 'real');
% 二关节数值示例：取消下一行注释并填写实际系数（应为非负值）。
% mat_Fv_value = [0.1; 0.2];
if ~(isa(mat_Fv_value, 'sym') || isnumeric(mat_Fv_value)) || ...
        ~isequal(size(mat_Fv_value), [N, 1])
    error('add_friction:InvalidCoefficient', '粘性摩擦系数必须是 N×1 符号或数值列向量。');
end

tau_before_friction = tau;
tau_viscous = mat_Fv_value .* dq_friction;
tau = simplify(tau_before_friction + tau_viscous);

%% 库伦摩擦（暂不启用）：取消以下代码注释即可追加到上述 tau
% mat_Fc_value = sym('Fc_', [N, 1], 'real'); % 每个关节的库伦摩擦系数，单位 N*m
% mat_Fc_value = [0.3; 0.4];               % 二关节数值示例，应填非负值
% tau_coulomb = mat_Fc_value .* sign(dq_friction);
% tau = simplify(tau + tau_coulomb);
% sign(0)=0；此模型不描述零速时的静摩擦。
