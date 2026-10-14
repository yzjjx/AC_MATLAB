% 从已提取的 M、C、G 与摩擦 F 构造 tau_a，对独立参数 pi_a 求雅可比得到 Y_a。
% 运行顺序（从工程根目录）：
%   addpath('MATLAB');
%   RNEA_SYM; Mat_out; add_friction; generate_Y_a;
% 沿用旧二连杆质心沿 x 轴的假设时：
%   RNEA_SYM; data_verify; Mat_out; add_friction; generate_Y_a;
% Mat_out 必须在 add_friction 前运行；不再需要 generate_Y_pi 或原始 Y。
%
% 原始参数质量/质心存在乘积，不能直接作为线性回归的独立求导变量。
% 本脚本提取物理参数单项式，将系数完全相同或成常数倍的项合并，
% 得到物理定义 pi_param，再用独立符号 pi_a 表示这些参数组合。
% 该过程不使用随机采样、SVD/QR；不保证任意模型都得到最小参数基。
%
% 主要输出：
%   pi_a、pi_param             独立参数符号及对应物理定义，p×1
%   M_a、C_a、G_a、F_a         改用 pi_a 表达的动力学矩阵/向量
%   tau_a_reference、tau_a     独立参考形式、代入跟踪误差后的控制补偿力矩
%   Y_a_reference、Y_a         分别由 jacobian(tau_a_*,pi_a) 直接得到
%   tau_a_reference_physical、tau_a_physical  原物理参数形式，供核对
%   e、de、v_r、a_r、s、Lambda               跟踪误差及参考量
% 保留前置脚本的 tau、M、C、G；F 使用实际 dq，不替换为 v_r。
% 关系：tau_a=Y_a*pi_a；subs(tau_a,pi_a,pi_param)=tau_a_physical。
% 数值控制使用 Y_a*pi_hat-Ks*s；pi_hat 顺序必须按本次 pi_param 排列。

%% 1. 输入与设置
required_a = {'N','mat_q','mat_dq','mat_ddq','M','C','G', ...
              'tau','tau_before_friction'};
for input_index_a = 1:numel(required_a)
    if ~exist(required_a{input_index_a},'var')
        error('generate_Y_a:MissingInput', ...
            '缺少 %s，请依次运行 RNEA_SYM、Mat_out、add_friction。', ...
            required_a{input_index_a});
    end
end
if ~isnumeric(N) || ~isscalar(N) || ~isfinite(N) || N < 1 || N ~= fix(N)
    error('generate_Y_a:InvalidJointCount','N 必须为正整数。');
end
if ~isa(mat_q,'sym') || ~isa(mat_dq,'sym') || ~isa(mat_ddq,'sym') || ...
        ~isequal(size(mat_q),[N,1]) || ~isequal(size(mat_dq),[N,1]) || ...
        ~isequal(size(mat_ddq),[N,1])
    error('generate_Y_a:InvalidMotion','q、dq、ddq 必须保留为 N×1 符号向量。');
end
motion_a = [mat_q;mat_dq;mat_ddq];
motion_names_a = cell(numel(motion_a),1);
for motion_index_a = 1:numel(motion_a)
    variables_a = symvar(motion_a(motion_index_a));
    if numel(variables_a) ~= 1 || ~isequal(variables_a,motion_a(motion_index_a))
        error('generate_Y_a:InvalidMotion','运动变量必须为独立符号。');
    end
    motion_names_a{motion_index_a} = char(motion_a(motion_index_a));
end
if numel(unique(motion_names_a)) ~= numel(motion_names_a)
    error('generate_Y_a:InvalidMotion','运动符号不能重复。');
end
if ~isa(M,'sym') || ~isa(C,'sym') || ~isa(G,'sym') || ...
        ~isa(tau,'sym') || ~isa(tau_before_friction,'sym') || ...
        ~isequal(size(M),[N,N]) || ~isequal(size(C),[N,N]) || ...
        ~isequal(size(G),[N,1]) || ~isequal(size(tau),[N,1]) || ...
        ~isequal(size(tau_before_friction),[N,1])
    error('generate_Y_a:InvalidDynamics','M、C、G 与 tau 的符号类型或维度不匹配。');
end
Lambda = 3*eye(N);                  % 常值误差增益；从这里修改
known_symbols_a = sym('g');         % 已知量；其余物理符号参与参数组合
display_full_matrix = true;

%% 2. 核对已提取矩阵，读取 add_friction 增加的实际速度力矩
if ~is_zero_a(tau_before_friction-M*mat_ddq-C*mat_dq-G)
    error('generate_Y_a:SourceMismatch', ...
        'M、C、G 与加入摩擦前的 tau 不一致，请重新运行同一模型的 Mat_out。');
end
if any(has(M,[mat_dq;mat_ddq]),'all') || ~is_zero_a(M-M.') || ...
        any(has(C,mat_ddq),'all') || any(has(G,[mat_dq;mat_ddq]),'all')
    error('generate_Y_a:InvalidDynamics','M、C、G 不符合当前刚体动力学结构。');
end
F = simplify(expand(tau-tau_before_friction));
if any(has(F,mat_ddq),'all') || ...
        ~is_zero_a(subs(F,mat_dq,-mat_dq)+F)
    error('generate_Y_a:InvalidFriction','F 必须为不含加速度的奇对称速度项。');
end

%% 3. 参考量及原物理参数形式的 tau_a
mat_qd = sym('qd_',[N,1],'real');
mat_dqd = sym('dqd_',[N,1],'real');
mat_ddqd = sym('ddqd_',[N,1],'real');
mat_vr = sym('vr_',[N,1],'real');
mat_ar = sym('ar_',[N,1],'real');
reference_symbols_a = [mat_qd;mat_dqd;mat_ddqd;mat_vr;mat_ar];
if any(has(motion_a,reference_symbols_a),'all') || ...
        any(has([M(:);C(:);G;F],reference_symbols_a),'all')
    error('generate_Y_a:SymbolCollision','前置模型不能使用 qd_、dqd_、ddqd_、vr_、ar_ 符号。');
end
e = mat_q-mat_qd;
de = mat_dq-mat_dqd;
v_r = mat_dqd-sym(Lambda)*e;
a_r = mat_ddqd-sym(Lambda)*de;
s = simplify(de+sym(Lambda)*e);
tau_a_reference_physical = simplify(M*mat_ar+C*mat_vr+G+F);
tau_a_physical = simplify(M*a_r+C*v_r+G+F);

%% 4. 更换参数：物理参数组合 -> 独立 pi_a
% 将 tau_a_reference 和动力学矩阵一起改写，保证可恢复全部 M、C、G、F。
% 首先处理力矩，使列顺序按其物理单项式的首次出现顺序确定。
physical_stack_a = [tau_a_reference_physical;M(:);C(:);G;F];
physical_symbols_a = setdiff(symvar(physical_stack_a).', ...
    [motion_a;reference_symbols_a;known_symbols_a(:)],'stable');
if isempty(physical_symbols_a)
    error('generate_Y_a:NoParameters','模型没有待估计参数，请保留物理参数为符号。');
end
[coefficient_groups_a,pi_param,parameter_terms_a] = ...
    group_parameters_a(physical_stack_a,physical_symbols_a);
parameter_count_a = numel(pi_param);
% 旧平面二连杆的七个物理定义完全匹配时，仅排序以兼容 AC_CON。
% 此步骤不构造 Y_a；参考矩阵仍在下一节通过雅可比求得。
if N == 2 && parameter_count_a == 7 && exist('mat_lc_value','var') && ...
        exist('mat_Izz_value','var') && exist('mat_Fv_value','var') && ...
        exist('mat_m','var') && exist('mat_L','var')
    preferred_definitions_a = [mat_Fv_value(1); ...
        mat_Izz_value(1)+mat_m(1)*mat_lc_value(1)^2+mat_m(2)*mat_L(1)^2; ...
        mat_m(2)*mat_L(1)*mat_lc_value(2); ...
        mat_m(1)*mat_lc_value(1)+mat_m(2)*mat_L(1); ...
        mat_Izz_value(2)+mat_m(2)*mat_lc_value(2)^2; ...
        mat_m(2)*mat_lc_value(2);mat_Fv_value(2)];
    order_a = zeros(7,1);
    for preferred_index_a = 1:7
        for definition_index_a = 1:7
            if is_zero_a(preferred_definitions_a(preferred_index_a)-pi_param(definition_index_a))
                order_a(preferred_index_a) = definition_index_a;
                break;
            end
        end
    end
    if all(order_a > 0) && numel(unique(order_a)) == 7
        coefficient_groups_a = coefficient_groups_a(:,order_a);
        pi_param = pi_param(order_a);
    end
end
pi_a = sym('pi_a_',[parameter_count_a,1],'real');
if any(has(physical_stack_a,pi_a),'all')
    error('generate_Y_a:SymbolCollision','前置模型使用了保留的 pi_a_ 参数符号。');
end
parameterized_stack_a = simplify(coefficient_groups_a*pi_a);
if ~is_zero_a(subs(parameterized_stack_a,pi_a,pi_param)-physical_stack_a)
    error('generate_Y_a:ParameterizationFailed','独立参数无法还原原物理模型。');
end
offset_a = N;
M_a = reshape(parameterized_stack_a(offset_a+(1:N*N)),N,N);
offset_a = offset_a+N*N;
C_a = reshape(parameterized_stack_a(offset_a+(1:N*N)),N,N);
offset_a = offset_a+N*N;
G_a = parameterized_stack_a(offset_a+(1:N));
offset_a = offset_a+N;
F_a = parameterized_stack_a(offset_a+(1:N));

%% 5. 直接构建 tau_a，并对独立参数求雅可比
tau_a_reference = simplify(M_a*mat_ar+C_a*mat_vr+G_a+F_a);
tau_a = simplify(M_a*a_r+C_a*v_r+G_a+F_a);
Y_a_reference = simplify(jacobian(tau_a_reference,pi_a));
Y_a = simplify(jacobian(tau_a,pi_a));
Y_friction_a = simplify(jacobian(F_a,pi_a));

%% 6. 验证线性、物理还原与实际运动量退化
reference_residual_a = simplify(expand(tau_a_reference-Y_a_reference*pi_a));
tracking_residual_a = simplify(expand(tau_a-Y_a*pi_a));
physical_residual_a = simplify(expand( ...
    subs(tau_a,pi_a,pi_param)-tau_a_physical));
reduction_residual_a = simplify(expand( ...
    subs(Y_a_reference,[mat_vr;mat_ar],[mat_dq;mat_ddq])*pi_param-tau));
substitution_residual_a = simplify(expand( ...
    subs(Y_a_reference,[mat_vr;mat_ar],[v_r;a_r])-Y_a));
if ~is_zero_a(reference_residual_a) || ~is_zero_a(tracking_residual_a) || ...
        ~is_zero_a(physical_residual_a) || ~is_zero_a(reduction_residual_a) || ...
        ~is_zero_a(substitution_residual_a) || ...
        any(has(Y_a,[pi_a;physical_symbols_a;mat_ddq]),'all')
    error('generate_Y_a:VerificationFailed','雅可比回归矩阵的符号验证失败。');
end

%% 7. 显示新参数定义与回归矩阵
fprintf('雅可比法 Y_a 已验证：%d 个关节，%d 个独立参数符号。\n',N,parameter_count_a);
disp('独立求导变量 pi_a ='); disp(pi_a);
disp('对应的物理参数组合 pi_param（列顺序以本次结果为准）='); disp(pi_param);
disp('Lambda ='); disp(Lambda);
disp('v_r ='); disp(v_r);
disp('a_r ='); disp(a_r);
disp('s ='); disp(s);
if display_full_matrix
    disp('M_a(pi_a) ='); disp(M_a);
    disp('C_a(q,dq,pi_a) ='); disp(C_a);
    disp('G_a(q,pi_a) ='); disp(G_a);
    disp('F_a(dq,pi_a) ='); disp(F_a);
    disp('tau_a(pi_a) ='); disp(tau_a);
    disp('Y_a_reference = jacobian(tau_a_reference,pi_a) ='); disp(Y_a_reference);
    disp('Y_a = jacobian(tau_a,pi_a) ='); disp(Y_a);
end
disp('验证通过：tau_a=Y_a*pi_a，代回 pi_param 后恢复物理力矩。');
% 数值导出（g 保留在 Y_a 中，导出前代入已知重力）：
%   Y_a_known = subs(Y_a,sym('g'),sym('9.81'));
%   Ya_fun = matlabFunction(Y_a_known, ...
%       'Vars',{mat_q,mat_dq,mat_qd,mat_dqd,mat_ddqd});
% 注意：原始 RNEA 模型保留 cy 等参数时，参数数目可能超过原七参数模型。
% 原 AC_CON 的硬编码七列只有在参数定义、列序与本次结果一致时才可使用。


function [groups,definitions,terms] = group_parameters_a(expressions,variables)
% 将多项式物理参数单项式收集为系数向量；仅合并完全相同/常数倍的系数。
count = numel(expressions);
terms = sym(zeros(0,1));
coefficients = sym(zeros(count,0));
names = {};
for row = 1:count
    try
        [values,row_terms] = coeffs(expand(expressions(row)),variables);
    catch exception
        error('generate_Y_a:NonpolynomialParameters', ...
            '物理参数须为多项式形式；第 %d 项无法提取：%s',row,exception.message);
    end
    for item = 1:numel(row_terms)
        if is_zero_a(values(item))
            continue;
        end
        if isequal(row_terms(item),sym(1))
            error('generate_Y_a:KnownTorqueOffset', ...
                '模型含与未知参数无关的非零项；请保持待估计动力学参数为符号。');
        end
        name = char(row_terms(item));
        [found,column] = ismember(name,names);
        if ~found
            column = numel(names)+1;
            names{column} = name; %#ok<AGROW>
            terms(column,1) = row_terms(item);
            coefficients(:,column) = sym(zeros(count,1));
        end
        coefficients(row,column) = coefficients(row,column)+values(item);
    end
end
groups = sym(zeros(count,0));
definitions = sym(zeros(0,1));
for column = 1:numel(terms)
    contribution = simplify(coefficients(:,column));
    if is_zero_a(contribution)
        continue;
    end
    merged = false;
    for group = 1:size(groups,2)
        pivot = find(~isAlways(groups(:,group)==0,'Unknown','false'),1);
        ratio = simplify(contribution(pivot)/groups(pivot,group));
        if isempty(symvar(ratio)) && is_zero_a(contribution-ratio*groups(:,group))
            definitions(group) = definitions(group)+ratio*terms(column);
            merged = true;
            break;
        end
    end
    if ~merged
        groups(:,end+1) = contribution; %#ok<AGROW>
        definitions(end+1,1) = terms(column); %#ok<AGROW>
    end
end
definitions = simplify(definitions);
end


function result = is_zero_a(expression)
result = isequal(simplify(expand(expression)),sym(zeros(size(expression))));
end
