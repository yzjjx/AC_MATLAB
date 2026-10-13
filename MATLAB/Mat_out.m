% 提取三矩阵M，C，G，输入为指定的tau，可以是简化之后，即通过data_verify.m处理之后的tau
% 标准动力学方程：tau = M*ddq + C*dq + G。
% 运行顺序：RNEA_SYM; Mat_out;
% 或：RNEA_SYM; data_verify; Mat_out;
% 可代入连杆参数数值，但应保留 mat_q、mat_dq、mat_ddq 对应的状态符号。
% 不清空工作区，不重新定义 N，也不修改原 tau。

required_vars = {'tau','N','mat_q','mat_dq','mat_ddq'};
for k = 1:numel(required_vars)
    if ~exist(required_vars{k}, 'var')
        error('Mat_out:MissingWorkspaceVariable', ...
            '工作区缺少变量 %s，请先运行 RNEA_SYM.m。', required_vars{k});
    end
end

%% 指定输入及要显示的矩阵
tau_input = tau; % 可改为工作区中的其他符号力矩列向量
output_matrices = {'M','C','G'}; % 例如只显示 M 时，改为 {'M'}

if ~isa(tau_input, 'sym') || ~isequal(size(tau_input), [N, 1]) || ...
        ~any(has(tau_input, mat_ddq), 'all')
    error('Mat_out:InvalidTau', ...
        '输入必须是保留加速度符号的 N×1 符号力矩列向量，不能用完全数值化的 tau 提取矩阵。');
end

%% 提取惯性矩阵：M(i,j) 是 tau_i 中 ddq_j 的系数
M = simplify(jacobian(tau_input, mat_ddq));

%% 提取重力向量：令速度、加速度均为零
G = simplify(subs(tau_input, [mat_dq;mat_ddq], sym(zeros(2*N, 1))));

%% 根据 Christoffel 公式构造标准科氏/离心矩阵
% C(i,j) = 1/2 * sum_k [dM(i,j)/dq_k + dM(i,k)/dq_j
%                      - dM(j,k)/dq_i] * dq_k
% C 的矩阵表示不唯一；这里采用与惯性矩阵 M 对应的标准构造。
M_derivative = sym(zeros(N, N, N));
for k = 1:N
    M_derivative(:,:,k) = diff(M, mat_q(k));
end
C = sym(zeros(N, N));
for i = 1:N
    for j = 1:N
        for k = 1:N
            C(i,j) = C(i,j) + sym(1)/2 * ...
                (M_derivative(i,j,k) + M_derivative(i,k,j) ...
                 - M_derivative(j,k,i)) * mat_dq(k);
        end
    end
end
C = simplify(C);

%% 检查分解是否能还原输入的力矩
tau_reconstructed = simplify(M*mat_ddq + C*mat_dq + G);
tau_residual = simplify(expand(tau_input - tau_reconstructed));
if ~isequal(tau_residual, sym(zeros(N, 1)))
    error('Mat_out:DecompositionMismatch', ...
        'M、C、G 未能还原输入 tau。请保留状态符号，并确认输入符合无摩擦刚体动力学方程。');
end

%% 按指定列表显示矩阵；M、C、G 均保存在当前工作区
for k = 1:numel(output_matrices)
    switch output_matrices{k}
        case 'M'
            disp('M ='); disp(M);
        case 'C'
            disp('C ='); disp(C);
        case 'G'
            disp('G ='); disp(G);
        otherwise
            error('Mat_out:UnknownMatrix', '输出名称只能为 M、C 或 G。');
    end
end
