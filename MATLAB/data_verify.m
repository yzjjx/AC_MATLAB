% RNEA_SYM.m 的延申脚本：读取当前工作区 tau，代入参数后仍输出 tau。
% 运行顺序：RNEA_SYM; data_verify;
% 不清空工作区，也不重新执行 RNEA。更改验证参数后应先重新运行 RNEA_SYM。

required_vars = {'tau','N','mat_c','mat_m','mat_L','mat_q','mat_dq', ...
                 'mat_ddq','I_xx','I_xy','I_xz','I_yy','I_yz','I_zz','g'};
for k = 1:numel(required_vars)
    if ~exist(required_vars{k}, 'var')
        error('data_verify:MissingWorkspaceVariable', ...
            '工作区缺少变量 %s，请先运行 RNEA_SYM.m。', required_vars{k});
    end
end
if ~isa(tau, 'sym') || ~isequal(size(tau), [N, 1])
    error('data_verify:InvalidTau', '工作区 tau 必须是 N×1 符号列向量。');
end

%% 验证参数：以下各项既可以填符号，也可以填数值
mat_lc_value = sym('lc_', [N, 1]); % 质心到关节原点的距离
mat_Izz_value = I_zz(:);          % 关于质心的 Z 轴惯量，尚未加入平行轴项
mat_m_value = mat_m;
mat_L_value = mat_L;
mat_q_value = mat_q;             % 当前 MDH 采用水平零位，不是论文竖直零位
mat_dq_value = mat_dq;
mat_ddq_value = mat_ddq;
g_value = g;

% 二连杆数值示例：需要时取消以下各行注释，覆盖上面的参数。
% 使用论文实物参数；q_1=0 对应论文 Q_1=pi/2，无摩擦力矩应为
% tau = [40.044982882; 0.083808] N*m。
% mat_lc_value = sym([0.091; 0.048]);
% mat_Izz_value = sym([1.266; 0.093]);
% mat_m_value = sym([23.902; 3.88]);
% mat_L_value = sym([0.45; 0.45]);
% mat_q_value = [sym(0); sym(pi)/2];
% mat_dq_value = sym([1; 2]);
% mat_ddq_value = sym([1; -1]);
% g_value = sym(9.81);

%==================================================================================无意义判断
parameter_vectors = {mat_lc_value, mat_Izz_value, mat_m_value, mat_L_value, ...
                     mat_q_value, mat_dq_value, mat_ddq_value};

for k = 1:numel(parameter_vectors)
    if ~(isa(parameter_vectors{k}, 'sym') || isnumeric(parameter_vectors{k})) || ...
            ~isequal(size(parameter_vectors{k}), [N, 1])
        error('data_verify:InvalidParameter', '验证参数向量必须均为 N×1。');
    end
end
if ~(isa(g_value, 'sym') || isnumeric(g_value)) || ~isscalar(g_value)
    error('data_verify:InvalidGravity', 'g_value 必须为符号或数值标量。');
end

%% 根据给定的平面验证假设建立目标质心和惯量
% 仅保留 Izz 是此处的验证设置，不代表真实三维惯量的其余分量必为零。
mat_c_value = [mat_lc_value.'; sym(zeros(2, N))];
mat_Ic_value = sym(zeros(3, 3, N));
for i = 1:N
    mat_Ic_value(:,:,i) = diag([sym(0), sym(0), mat_Izz_value(i)]);
    % 生成 I_c1、I_c2、... 和 c_1、c_2、...，供直接查看
    eval(sprintf('I_c%d = mat_Ic_value(:,:,%d);', i, i));
    eval(sprintf('c_%d = mat_c_value(:,%d);', i, i));
end

%% 对工作区已有 tau 做符号代入；重新赋值 mat_Ic 不会改变已生成的 tau
tau_before_verify = tau;
old_parameters = [I_xx(:); I_xy(:); I_xz(:); I_yy(:); I_yz(:); I_zz(:); ...
                  mat_c(:); mat_m; mat_L; mat_q; mat_dq; mat_ddq; g];
new_parameters = [sym(zeros(5*N, 1)); mat_Izz_value; mat_c_value(:); ...
                  mat_m_value; mat_L_value; mat_q_value; ...
                  mat_dq_value; mat_ddq_value; g_value];
tau = simplify(subs(tau_before_verify, old_parameters, new_parameters));

% 命令行tau(1)为tau_1，以此类推