% 复用 generate_Y_pi 已得到的 pi_param，构造 tau_a 并求雅可比回归矩阵 Y_a。
% 从工程根目录运行：addpath('MATLAB');
% RNEA_SYM; Mat_out; add_friction; generate_Y_pi; generate_Y_a;
% 质心沿 x 轴时：RNEA_SYM; data_verify; Mat_out; add_friction; generate_Y_pi; generate_Y_a;
% Mat_out 在加摩擦前运行；generate_Y_pi 在加摩擦后运行，使 pi_param 包含摩擦参数。
%
% 输入：已有的 N、q/dq/ddq、M/C/G、tau、tau_before_friction，及 generate_Y_pi 的
% pi_param、parameter_symbols、term_names、base_columns（均须来自同一模型）。
% 本脚本不重新求解、合并或排序物理参数；保持 pi_param 及其列顺序不变。
% pi_a 只是已有参数组合的独立符号表示，用于 jacobian 求导。
% 输出：M_a/C_a/G_a/F_a、tau_a、Y_a；tau_a=Y_a*pi_a。
% 代入 pi_a=pi_param 后恢复物理模型；F 始终使用实际速度 dq。

%% 1. 设置增益，构造参考运动量
Lambda = 3*eye(N);                          % 常值误差增益 Lambda=3I。
F = simplify(tau-tau_before_friction);       % 读取加入摩擦前后的力矩差。
mat_qd = sym('qd_',[N,1],'real');            % 期望位置符号。
mat_dqd = sym('dqd_',[N,1],'real');           % 期望速度符号。
mat_ddqd = sym('ddqd_',[N,1],'real');         % 期望加速度符号。
mat_vr = sym('vr_',[N,1],'real');            % 独立参考速度，用于显示一般形式。
mat_ar = sym('ar_',[N,1],'real');            % 独立参考加速度。
e = mat_q-mat_qd;                           % 位置误差 e=q-qd。
de = mat_dq-mat_dqd;                        % 速度误差 de=dq-dqd。
v_r = mat_dqd-sym(Lambda)*e;                % 参考速度 vr=dqd-Lambda*e。
a_r = mat_ddqd-sym(Lambda)*de;              % 参考加速度 ar=ddqd-Lambda*de。
s = simplify(de+sym(Lambda)*e);             % 滤波误差 s=de+Lambda*e=dq-vr。

%% 2. 用已有 pi_param 的独立符号表示 M、C、G、F
parameter_count_a = numel(pi_param);         % 直接读取前置脚本得到的参数数量。
pi_a = sym('pi_a_',[parameter_count_a,1],'real'); % 仅创建独立符号，不计算新的物理参数组合。
physical_stack_a = [M(:);C(:);G;F];         % 将四个矩阵/向量按列堆叠，统一更换参数表示。
base_term_names_a = term_names(base_columns); % 前置脚本选定的主元单项式，顺序与 pi_param 一致。
coefficient_stack_a = sym(zeros(numel(physical_stack_a),parameter_count_a)); % 各参数的系数矩阵。

% generate_Y_pi 的 parameter_map(:,base_columns)=I，因此各主元单项式的系数
% 就是相应 pi_a 的系数。这里只提取已有基的系数，不再判断秩或合并参数。
for row_a = 1:numel(physical_stack_a)        % 每次处理一个矩阵元素。
    [values_a,terms_a] = coeffs(expand(physical_stack_a(row_a)),parameter_symbols); % 沿用前置物理符号。
    for item_a = 1:numel(terms_a)           % 遍历当前元素的物理参数单项式。
        [found_a,column_a] = ismember(char(terms_a(item_a)),base_term_names_a); % 查找已有主元列号。
        if found_a                         % 只保留主元单项式的系数，其余项由原参数组合统一表示。
            coefficient_stack_a(row_a,column_a) = ...
                coefficient_stack_a(row_a,column_a)+values_a(item_a); % 累加相应参数的系数。
        end
    end
end
parameterized_stack_a = simplify(coefficient_stack_a*pi_a); % 用独立符号替代已有参数组合。

% 按上面的堆叠顺序恢复矩阵；reshape 与 (:) 的按列展开互为逆操作。
offset_a = 0;                              % 从长向量开头读取。
M_a = reshape(parameterized_stack_a(offset_a+(1:N*N)),N,N); % 前 N^2 项：惯性矩阵。
offset_a = offset_a+N*N;                    % 游标移到 M 之后。
C_a = reshape(parameterized_stack_a(offset_a+(1:N*N)),N,N); % 再取 N^2 项：科氏/离心矩阵。
offset_a = offset_a+N*N;                    % 游标移到 C 之后。
G_a = parameterized_stack_a(offset_a+(1:N)); % 接下来 N 项：重力向量。
offset_a = offset_a+N;                      % 游标移到 G 之后。
F_a = parameterized_stack_a(offset_a+(1:N)); % 最后 N 项：摩擦向量。

%% 3. 从 M_a、C_a、G_a、F_a 构造 tau_a，对已有参数的独立符号求雅可比
% C_a 仍在实际 q、dq 处计算，只把右乘向量换成 vr；F_a 保留实际 dq。
tau_a_reference_physical = simplify(M*mat_ar+C*mat_vr+G+F); % 原物理参数形式，供核对。
tau_a_physical = simplify(M*a_r+C*v_r+G+F);  % 代入误差反馈后的原物理力矩。
tau_a_reference = simplify(M_a*mat_ar+C_a*mat_vr+G_a+F_a); % 独立参考量形式。
tau_a = simplify(M_a*a_r+C_a*v_r+G_a+F_a);   % 独立参数形式的控制补偿力矩。
Y_a_reference = simplify(jacobian(tau_a_reference,pi_a)); % 对已有参数的独立符号求导。
Y_a = simplify(jacobian(tau_a,pi_a));        % 所需回归矩阵，第 j 列为对 pi_a(j) 的偏导。
Y_friction_a = simplify(jacobian(F_a,pi_a)); % 单独查看摩擦部分的回归矩阵。

%% 4. 验证已有参数基能够恢复各矩阵和原力矩
parameterization_residual_a = simplify(expand( ...
    subs(parameterized_stack_a,pi_a,pi_param)-physical_stack_a)); % 换回已有物理定义，检查全部矩阵。
reference_residual_a = simplify(expand(tau_a_reference-Y_a_reference*pi_a)); % 检查参考形式的线性关系。
tracking_residual_a = simplify(expand(tau_a-Y_a*pi_a)); % 检查 tau_a=Y_a*pi_a。
reduction_residual_a = simplify(expand( ...
    subs(Y_a_reference,[mat_vr;mat_ar],[mat_dq;mat_ddq])*pi_param-tau)); % 实际参考量下应恢复 tau。
verification_residual_a = [parameterization_residual_a;reference_residual_a; ...
    tracking_residual_a;reduction_residual_a]; % 汇总已化简的残差。
if ~isequal(verification_residual_a,sym(zeros(size(verification_residual_a)))) || ...
        any(has(Y_a,[pi_a;parameter_symbols;mat_ddq]),'all') % Y_a 不应包含待估计参数或实际加速度。
    error('generate_Y_a:VerificationFailed', ...
        '已有参数基无法还原当前模型，请对同一模型依次运行 Mat_out、add_friction、generate_Y_pi。');
end

%% 5. 显示结果（所有变量同时保留在工作区）
fprintf('Y_a 已验证：%d 个关节，复用 %d 个已有参数组合。\n',N,parameter_count_a); % 显示模型规模。
disp('已有物理参数 pi_param（保持 generate_Y_pi 的顺序）='); disp(pi_param);
disp('对应的独立符号 pi_a ='); disp(pi_a);
disp('Lambda ='); disp(Lambda);
disp('v_r ='); disp(v_r);
disp('a_r ='); disp(a_r);
disp('s ='); disp(s);
disp('M_a ='); disp(M_a);
disp('C_a ='); disp(C_a);
disp('G_a ='); disp(G_a);
disp('F_a ='); disp(F_a);
disp('tau_a ='); disp(tau_a);
disp('Y_a_reference ='); disp(Y_a_reference);
disp('Y_a ='); disp(Y_a);
disp('验证通过：tau_a=Y_a*pi_a，代入已有 pi_param 后恢复物理模型。');

% 数值控制：Y_a*pi_hat-Ks*s；pi_hat 必须沿用本次 pi_param 的顺序。
% AC_CON.m 内的七列顺序为固定顺序，直接使用前需按物理定义核对列的对应关系。
% 导出示例：Y_a_known = subs(Y_a,sym('g'),sym('9.81'));
% Ya_fun = matlabFunction(Y_a_known,'Vars',{mat_q,mat_dq,mat_qd,mat_dqd,mat_ddqd});
