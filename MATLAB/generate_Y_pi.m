% 从工作区的符号 tau 提取回归矩阵与固定参数组合：tau = Y*pi_param。
% 流程：随机运动状态 -> 采样符号力矩 -> pi_param -> 反推符号 Y。
% 方法：SVD 判断独立参数数量，QR 选择独立列和独立采样行。
% 运行顺序：RNEA_SYM; generate_Y_pi;
% 也可以先运行 data_verify，但 q、dq、ddq 和待辨识参数必须保留为符号。
% 本脚本不修改 tau；pi_param 表示参数向量，避免覆盖 MATLAB 的常数 pi。

% 无意义判断
% =============================================================================================
if ~exist('tau', 'var') || ~isa(tau, 'sym') || ~iscolumn(tau)
    error('generate_Y_pi:InvalidTau', '请先在工作区生成符号列向量 tau。');
end
if ~exist('mat_q', 'var') || ~exist('mat_dq', 'var') || ~exist('mat_ddq', 'var')
    error('generate_Y_pi:MissingMotion', '工作区需要 mat_q、mat_dq、mat_ddq。');
end
tau_for_Y = tau;
joint_count = numel(tau_for_Y);
motion_symbols = [mat_q(:); mat_dq(:); mat_ddq(:)];

%% 设置：只随机代入运动状态，不随机代入待辨识参数
% 默认 g 是已知常量：采样时取 9.81，最终 Y 中仍保留符号 g。
% 长度默认作为未知参数进入 pi_param，与示例中的 m2*l1 一致。
known_symbols = sym('g');
known_values = 9.81;
% 如果连杆长度已知，可改为（数值按实际模型填写）：
% known_symbols = [g; mat_L];
% known_values = [9.81; 0.45; 0.45];

random_seed = 0;
rank_tol = 1e-10;                 % 列归一化后的相对秩阈值
parameter_symbols = setdiff(symvar(tau_for_Y).', ...
                            [motion_symbols; known_symbols], 'stable');
if isempty(parameter_symbols)
    error('generate_Y_pi:NoParameters', 'tau 中没有保留符号形式的待辨识参数。');
end

%% 1. 先随机代入 q、dq、ddq，得到仍含未知物理参数的 tau_samples
% 采样次数取 20 与参数符号数量四倍中的较大值；符号包括 m、L、lc、Izz 等。
n_samples = max(20, 4*numel(parameter_symbols)); % 可增加采样数量

% sample_stream为随机数生成器，随机数生成算法为mt19937ar，随机数种子为random_seed
sample_stream = RandStream('mt19937ar', 'Seed', random_seed);

% 令 q = 2*atan(有理数)，使 sin(q)、cos(q) 可精确化为有理数。
% 这样既能随机采样，也避免用浮点数拼出近似的参数组合。
half_angle_tan = sym(randi(sample_stream, [-100, 100], ...
                           joint_count, n_samples))/20;
q_samples_exact = 2*atan(half_angle_tan);
% 生成随机速度
dq_samples_exact = sym(randi(sample_stream, [-6, 6], ...
                            joint_count, n_samples))/2;
% 生成随机加速度
ddq_samples_exact = sym(randi(sample_stream, [-10, 10], ...
                             joint_count, n_samples))/2;
% 生成随机运动状态矩阵，每列为一次实验的 q、dq、ddq 
motion_samples_exact = [q_samples_exact; dq_samples_exact; ddq_samples_exact];

% 转换为double，可以绘图显示这个q和dq和ddq这些随机数
q_samples = double(q_samples_exact);
dq_samples = double(dq_samples_exact);
ddq_samples = double(ddq_samples_exact);

% 创建符号矩阵，用来保存每一次计算的力矩 tau，仍然含有未知参数
tau_samples = sym(zeros(joint_count, n_samples));

for sample_idx = 1:n_samples
    % 每列对应一次实验；这里直接对工作区的原始 tau 进行代入。
    tau_samples(:, sample_idx) = simplify(subs(tau_for_Y, motion_symbols, ...
                                              motion_samples_exact(:, sample_idx)));
end

%% 2. 从采样力矩提取固定参数组合 pi_param（此时还没有生成符号 Y）
% theta_full 存放采样力矩中出现的物理参数单项式，例如 m1、m2*l1。
% W_symbolic*theta_full 等于 tau_samples 按列堆叠后的力矩向量。

% 存储动力学方程出现的全部参数单项式，theta_full
theta_full = sym(zeros(0, 1));

% N*S,N为关节数，S为随机输入次数，存储每一个参数单项式对应的系数
W_symbolic = sym(zeros(joint_count*n_samples, 0));

% 记录已经出现过哪些参数单项式，防止重复添加
term_names = {};

% 开始循环每次采样n_samples，每个关节joint_count
for sample_idx = 1:n_samples
    for joint_idx = 1:joint_count
        row_idx = (sample_idx-1)*joint_count + joint_idx;
        % expand用于展开所有计算，tau_samples为2行28列的数据
        % coeffs用来分解成系数和单项式，[C, T] = coeffs(tau, [m1,m2,lc2,Izz1]);
        % 例如2m_1，就会被分解成C为2，T为m_1
        [term_coefficients, parameter_terms] = ...
            coeffs(expand(tau_samples(joint_idx, sample_idx)), parameter_symbols);
        
        for term_idx = 1:numel(parameter_terms)
            % 转换为char字符串，因为同一个参数可能出现在不同关节和采样中
            % ismember用来确定 A 的哪些行同时也在 B 中及其在 B 中的相应位置
            % 例如 m_1 在 28 次采样中可能出现很多次，但 theta_full 中只需要记录一次
            term_name = char(parameter_terms(term_idx));
            [found, column_idx] = ismember(term_name, term_names);

            % 当参数第一次出现，代码就为这个参数创建新的一列
            if ~found
                column_idx = numel(term_names) + 1;
                term_names{column_idx} = term_name;
                theta_full(column_idx, 1) = parameter_terms(term_idx);
                W_symbolic(:, column_idx) = sym(zeros(joint_count*n_samples, 1));
            end

            % 然后将对应的系数放到这个矩阵里面
            W_symbolic(row_idx, column_idx) = ...
                W_symbolic(row_idx, column_idx) + term_coefficients(term_idx);
        end
    end
end
parameter_count = numel(theta_full);

% =================最小参数集的提取==========================
% SVD 用来解释有多少个独立参数组合，QR 用来选择哪几列作为基，而符号线性方程求解用来表示这些参数如何组合
% SVD 判断数值秩，QR 选择主元；参数关系随后由精确符号运算得到。

% 首先代入已知的数值参数，如重力等，生成数值矩阵 W_scaled，用于 SVD 和 QR
W = double(subs(W_symbolic, known_symbols, known_values));
% 列归一化：每列除以其二范数，避免量纲差异导致的数值秩判断错误
% 假如W = [1000 0.001;2000 0.002],避免第一列的数值过大导致第二列被忽略，归一化后W_scaled = [1 1;2 2]
column_scale = vecnorm(W, 2, 1);
column_scale(column_scale == 0) = 1;
% 列归一化，通过这种归一化的方式，每列的二范数为1
W_scaled = W ./ column_scale;

% SVD：返回按降序排列的奇异值，超过相对阈值的数量就是数值秩。
% 数值秩 n_base 决定需要保留多少个独立参数组合。
singular_values = svd(W_scaled);
% 如果第三列可以由前面两列表示，则奇异值第三个为0，n_base=2

if isempty(singular_values) || singular_values(1) == 0
    svd_threshold = 0;
    n_base = 0;
else
    % 设置阈值进行判断，如果奇异值小于阈值，则认为是线性相关的，n_base为独立参数组合的数量
    svd_threshold = rank_tol*singular_values(1);
    n_base = sum(singular_values > svd_threshold);
end

% 带列置换的 QR：按独立性选择原矩阵的列，保留前 n_base 个列号。
% 假设SVD已经找出矩阵秩为5，那么在原来的N列里面，选择哪5列作为独立参数组合
% QR分解的列置换可以帮助选择这些独立的列
[~, R_sample, column_order] = qr(W_scaled, 'vector');
base_columns = sort(column_order(1:n_base));

if n_base == 0
    parameter_map = sym(zeros(0, parameter_count));
else
    % 第二次 QR：转置后选列，等价于在原矩阵中选独立采样行。
    % 选出 n_base 行，使主元列对应的子矩阵成为可逆方阵。
    [~, ~, row_order] = qr(W_scaled(:, base_columns).', 'vector');
    exact_rows = row_order(1:n_base);
    A_exact = W_symbolic(exact_rows, :);
    % pi_param = parameter_map*theta_full；主元列构成单位矩阵。
    % g 和已知几何符号在这里不被数值替换。
    parameter_map = simplify(A_exact(:, base_columns) \ A_exact);
end
pi_param = simplify(parameter_map*theta_full);

%% 3. 再回到原始符号 tau，通过系数匹配反推符号回归矩阵 Y
% tau = Y_full*theta_full = Y*parameter_map*theta_full。
% 因 parameter_map(:,base_columns)=eye(n_base)，
% 对这些主元单项式匹配系数，即得到 Y = Y_full(:,base_columns)。
Y_full = sym(zeros(joint_count, parameter_count));
for joint_idx = 1:joint_count
    [term_coefficients, parameter_terms] = ...
        coeffs(expand(tau_for_Y(joint_idx)), parameter_symbols);
    for term_idx = 1:numel(parameter_terms)
        [found, column_idx] = ismember(char(parameter_terms(term_idx)), term_names);
        if ~found && ~isequal(term_coefficients(term_idx), sym(0))
            error('generate_Y_pi:InsufficientSamples', ...
                '采样遗漏了原始 tau 中的参数项，请增加 n_samples。');
        elseif found
            Y_full(joint_idx, column_idx) = ...
                Y_full(joint_idx, column_idx) + term_coefficients(term_idx);
        end
    end
end
Y_full = simplify(Y_full);
Y = Y_full(:, base_columns);

%% 4. 符号验证：参数组合固定，且对任意运动状态满足原动力学
regressor_residual = simplify(expand(Y_full - Y*parameter_map));
tau_Y_pi = simplify(Y*pi_param);
tau_residual = simplify(expand(tau_for_Y - tau_Y_pi));
if any(has(pi_param, motion_symbols), 'all') || ...
        ~isequal(regressor_residual, sym(zeros(size(Y_full)))) || ...
        ~isequal(tau_residual, sym(zeros(size(tau_for_Y))))
    error('generate_Y_pi:VerificationFailed', ...
        '符号验证未通过，请增加 n_samples 或减小 rank_tol 后重新运行。');
end

fprintf('采样次数：%d；参数单项式：%d；独立参数组合：%d。\n', ...
        n_samples, parameter_count, n_base);
disp('固定参数向量 pi_param =');
disp(pi_param);
disp('反推得到的符号回归矩阵 Y =');
disp(Y);
disp('前三次随机实验的符号力矩（每列为一次实验）=');
disp(tau_samples(:, 1:min(3, n_samples)));
disp('验证通过：tau - Y*pi_param = 0。');
