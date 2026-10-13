% 从现有 Y 和 pi_param 生成普通 AC 的控制回归矩阵 Y_a。
% 适用：固定基座、全旋转关节的刚体串联机器人，关节数不写死。
% 前置模型必须正确；本脚本不修改 MDH、不重新建模、不重新提取参数基。
%
% 示例运行顺序（从本工程根目录运行）：
%   addpath(fullfile(pwd, 'MATLAB'));
%   RNEA_SYM; data_verify; generate_Y_pi; generate_Y_a;
% 若加入摩擦：
%   RNEA_SYM; data_verify; add_friction; generate_Y_pi; generate_Y_a;
% 不需要先运行 Mat_out；本脚本按 Y 的每一列提取动力学贡献。
%
% 输入：Y、pi_param、mat_q、mat_dq、mat_ddq（均保留符号）。
% 输出：
%   Lambda                         常值对称正定矩阵
%   mat_qd/mat_dqd/mat_ddqd         期望位置/速度/加速度符号
%   e、de、v_r、a_r、s              图中误差和参考量
%   mat_vr、mat_ar                  独立参考速度/加速度符号
%   Y_a_reference                  Y_a(q,dq,vr,ar)，N×p
%   Y_a                            代入 v_r、a_r 后的 Y_a，N×p
%   M_a、C_a、G_a、F_a              同一参数基重构的动力学各项
%   tau_a_reference、tau_a          控制补偿表达式，不包含 -Ks*s
%
% Y_a*pi_param = M_a*a_r + C_a*v_r + G_a + F_a。
% F_a 为可选黏性/库伦等奇对称速度项，使用实际 dq，不替换为 v_r。
% 在线控制：tau_cmd = Y_a_numeric*pi_hat - Ks*s_numeric。
% 本脚本只生成符号矩阵，不读取轨迹 CSV，也不发送设备指令。

%% 1. 显式设置 Lambda（从这里修改增益）
if ~exist('mat_q', 'var') || ~isa(mat_q, 'sym') || ...
        ~iscolumn(mat_q) || isempty(mat_q)
    error('generate_Y_a:InvalidPosition', ...
        '请先生成非空符号列向量 mat_q，再运行 generate_Y_a。');
end

% 对任意 N：Lambda = diag(lambda_1,...,lambda_N)，单位 s^(-1)。
% 默认每个关节为 3；二关节时就是 [3,0; 0,3]。
% 不同关节不同增益可改为 Lambda = diag([3, 4, ...]);
Lambda = 3 * eye(numel(mat_q));

% true：允许当前 add_friction.m 的黏性/库伦摩擦等奇对称速度项。
% false：要求 Y 只包含无摩擦刚体动力学。
include_friction = true;
display_full_matrix = true;

%% 2. 检查前置模型，保持现有参数基与列顺序
required_vars_a = {'Y', 'pi_param', 'mat_dq', 'mat_ddq'};
for check_idx_a = 1:numel(required_vars_a)
    if ~exist(required_vars_a{check_idx_a}, 'var')
        error('generate_Y_a:MissingInput', ...
            '工作区缺少 %s，请先运行建模和 generate_Y_pi。', ...
            required_vars_a{check_idx_a});
    end
end
joint_count_a = numel(mat_q);
if exist('N', 'var') && (~isnumeric(N) || ~isscalar(N) || N ~= joint_count_a)
    error('generate_Y_a:JointCountMismatch', 'N 必须与 mat_q 的长度一致。');
end
N = joint_count_a;
if ~isa(mat_dq, 'sym') || ~isa(mat_ddq, 'sym') || ...
        ~isequal(size(mat_dq), [N,1]) || ~isequal(size(mat_ddq), [N,1])
    error('generate_Y_a:InvalidMotion', ...
        'mat_q、mat_dq、mat_ddq 必须是同长度的符号列向量。');
end
if ~isa(Y, 'sym') || size(Y,1) ~= N || isempty(Y) || ...
        ~isa(pi_param, 'sym') || ~iscolumn(pi_param) || ...
        numel(pi_param) ~= size(Y,2)
    error('generate_Y_a:InvalidRegressor', ...
        'Y 必须是 N×p 符号矩阵，pi_param 必须是配套的 p×1 符号列向量。');
end
motion_symbols_a = [mat_q; mat_dq; mat_ddq];
motion_names_a = cell(numel(motion_symbols_a),1);
for check_idx_a = 1:numel(motion_symbols_a)
    scalar_vars_a = symvar(motion_symbols_a(check_idx_a));
    if numel(scalar_vars_a) ~= 1 || ...
            ~isequal(motion_symbols_a(check_idx_a), scalar_vars_a)
        error('generate_Y_a:InvalidMotionSymbol', ...
            '运动变量必须为独立符号，不能提前代入数值或表达式。');
    end
    motion_names_a{check_idx_a} = char(motion_symbols_a(check_idx_a));
end
if numel(unique(motion_names_a)) ~= numel(motion_names_a) || ...
        any(has(pi_param, motion_symbols_a), 'all')
    error('generate_Y_a:InvalidParameterBasis', ...
        '运动符号必须互不重复，pi_param 必须不依赖运动状态。');
end
if ~isnumeric(Lambda) || ~isreal(Lambda) || ...
        ~isequal(size(Lambda), [N,N]) || any(~isfinite(Lambda), 'all')
    error('generate_Y_a:InvalidLambda', 'Lambda 必须是有限实数 N×N 矩阵。');
end
if ~isequal(Lambda, Lambda.')
    error('generate_Y_a:InvalidLambda', 'Lambda 必须为对称矩阵。');
end
[~, lambda_chol_flag_a] = chol(Lambda);
if lambda_chol_flag_a ~= 0
    error('generate_Y_a:InvalidLambda', 'Lambda 必须为正定矩阵。');
end

% 检查工作区 tau 是否与本次 Y 对应，避免混用不同模型的旧变量。
tau_from_Y_a = simplify(Y*pi_param);
if exist('tau', 'var')
    if ~isa(tau, 'sym') || ~isequal(size(tau), [N,1]) || ...
            ~isequal(simplify(expand(tau-tau_from_Y_a)), sym(zeros(N,1)))
        error('generate_Y_a:SourceMismatch', ...
            'tau 与 Y*pi_param 不一致，请重新生成同一模型的 Y 和 pi_param。');
    end
end

%% 3. 显式呈现图中的公式：期望量、误差、参考量
mat_qd = sym('qd_', [N,1], 'real');
mat_dqd = sym('dqd_', [N,1], 'real');
mat_ddqd = sym('ddqd_', [N,1], 'real');
if any(has(motion_symbols_a, [mat_qd;mat_dqd;mat_ddqd]), 'all')
    error('generate_Y_a:ReferenceSymbolCollision', ...
        '实际运动变量不能与期望运动符号 qd_、dqd_、ddqd_ 重名。');
end

e = mat_q - mat_qd;                   % e = q - q_d
de = mat_dq - mat_dqd;                % de = dq - dq_d
v_r = mat_dqd - sym(Lambda)*e;        % v_r = dq_d - Lambda*e
a_r = mat_ddqd - sym(Lambda)*de;      % a_r = ddq_d - Lambda*de
s = simplify(mat_dq - v_r);          % s = dq - v_r = de + Lambda*e

% 独立参考符号用于先验证一般恒等式，随后再代入上述误差公式。
mat_vr = sym('vr_', [N,1], 'real');
mat_ar = sym('ar_', [N,1], 'real');
if any(has(motion_symbols_a, [mat_vr;mat_ar]), 'all')
    error('generate_Y_a:ReferenceSymbolCollision', ...
        '实际运动变量不能与独立参考符号 vr_、ar_ 重名。');
end

%% 4. 按原 Y 逐列构造：Y_a_reference(:,j) = Mj*ar + Cj*vr + Gj + Fj
parameter_count_a = size(Y,2);
Y_a_reference = sym(zeros(N, parameter_count_a));
Y_friction_a = sym(zeros(N, parameter_count_a));
M_a = sym(zeros(N));
C_a = sym(zeros(N));
G_a = sym(zeros(N,1));
F_a = sym(zeros(N,1));

for column_idx_a = 1:parameter_count_a
    y_column_a = Y(:,column_idx_a);
    M_column_a = simplify(jacobian(y_column_a, mat_ddq));
    if any(has(M_column_a, [mat_dq;mat_ddq]), 'all') || ...
            ~isequal(simplify(M_column_a-M_column_a.'), sym(zeros(N)))
        error('generate_Y_a:InvalidInertia', ...
            'Y 第 %d 列的惯性贡献须对称且只依赖位置和固定已知量。', column_idx_a);
    end
    G_column_a = simplify(subs(y_column_a, [mat_dq;mat_ddq], sym(zeros(2*N,1))));

    % 与 Mat_out 一致：Christoffel 构造。C 内部保持实际 mat_dq。
    M_derivatives_a = sym(zeros(N,N,N));
    for k_a = 1:N
        M_derivatives_a(:,:,k_a) = diff(M_column_a, mat_q(k_a));
    end
    C_column_a = sym(zeros(N));
    for i_a = 1:N
        for j_a = 1:N
            for k_a = 1:N
                C_column_a(i_a,j_a) = C_column_a(i_a,j_a) + sym(1)/2 * ...
                    (M_derivatives_a(i_a,j_a,k_a) ...
                    +M_derivatives_a(i_a,k_a,j_a) ...
                    -M_derivatives_a(j_a,k_a,i_a))*mat_dq(k_a);
            end
        end
    end
    C_column_a = simplify(C_column_a);

    % 例如黏性摩擦列只留下 dq_i；库伦摩擦列留下 sign(dq_i)。
    F_column_a = simplify(expand(y_column_a - M_column_a*mat_ddq ...
                                - C_column_a*mat_dq - G_column_a));
    if ~isequal(F_column_a, sym(zeros(N,1)))
        if ~include_friction
            error('generate_Y_a:UnexpectedVelocityTerm', ...
                'Y 第 %d 列含附加速度项，请检查模型或启用 include_friction。', column_idx_a);
        end
        odd_residual_a = simplify(expand( ...
            subs(F_column_a, mat_dq, -mat_dq) + F_column_a));
        if any(has(F_column_a, mat_ddq), 'all') || ...
                ~isequal(odd_residual_a, sym(zeros(N,1)))
            error('generate_Y_a:UnsupportedDynamics', ...
                ['Y 第 %d 列不符合刚体动力学加奇对称摩擦项；' ...
                 '不能将未匹配的惯性/科氏项当作摩擦。'], column_idx_a);
        end
    end

    Y_a_reference(:,column_idx_a) = simplify( ...
        M_column_a*mat_ar + C_column_a*mat_vr + G_column_a + F_column_a);
    Y_friction_a(:,column_idx_a) = F_column_a;
    M_a = M_a + M_column_a*pi_param(column_idx_a);
    C_a = C_a + C_column_a*pi_param(column_idx_a);
    G_a = G_a + G_column_a*pi_param(column_idx_a);
    F_a = F_a + F_column_a*pi_param(column_idx_a);
end
M_a = simplify(M_a);
C_a = simplify(C_a);
G_a = simplify(G_a);
F_a = simplify(F_a);

%% 5. 验证一般参考回归关系，再代入图中的 v_r 和 a_r
tau_a_reference = simplify(M_a*mat_ar + C_a*mat_vr + G_a + F_a);
Y_a = simplify(subs(Y_a_reference, [mat_vr;mat_ar], [v_r;a_r]));
tau_a = simplify(M_a*a_r + C_a*v_r + G_a + F_a);

source_residual_a = simplify(expand( ...
    tau_from_Y_a - M_a*mat_ddq - C_a*mat_dq - G_a - F_a));
reference_residual_a = simplify(expand(tau_a_reference - Y_a_reference*pi_param));
tracking_residual_a = simplify(expand(tau_a - Y_a*pi_param));
reduction_residual_a = simplify(expand( ...
    subs(Y_a_reference, [mat_vr;mat_ar], [mat_dq;mat_ddq]) - Y));
if ~isequal(source_residual_a, sym(zeros(N,1))) || ...
        ~isequal(reference_residual_a, sym(zeros(N,1))) || ...
        ~isequal(tracking_residual_a, sym(zeros(N,1))) || ...
        ~isequal(reduction_residual_a, sym(zeros(size(Y)))) || ...
        any(has(Y_a, mat_ddq), 'all')
    error('generate_Y_a:VerificationFailed', '控制回归矩阵的符号恒等式验证未通过。');
end

%% 6. 显示结果；保持输入 tau、Y、pi_param、M、C、G 不变
fprintf('Y_a 已生成并验证：%d 个旋转关节，%d 个参数组合。\n', N, parameter_count_a);
disp('Lambda ='); disp(Lambda);
disp('e = q - q_d ='); disp(e);
disp('de = dq - dq_d ='); disp(de);
disp('v_r = dq_d - Lambda*e ='); disp(v_r);
disp('a_r = ddq_d - Lambda*de ='); disp(a_r);
disp('s = dq - v_r = de + Lambda*e ='); disp(s);
if display_full_matrix
    disp('Y_a_reference(q,dq,vr,ar) ='); disp(Y_a_reference);
    disp('Y_a(q,dq,qd,dqd,ddqd)，已代入 Lambda 和误差公式 ='); disp(Y_a);
end
disp('验证通过：Y_a*pi_param = M_a*a_r + C_a*v_r + G_a + F_a。');
% 不需要额外一阶低通滤波器来生成本矩阵；s 是误差组合，不是测量滤波。
% 已知几何量/重力数值化后，可用 matlabFunction 导出数值函数：
%   matlabFunction(Y_a_known,'Vars',{mat_q,mat_dq,mat_qd,mat_dqd,mat_ddqd});
% 或导出 Y_a_reference_known，输入 {mat_q,mat_dq,mat_vr,mat_ar}。
