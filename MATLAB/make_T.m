function T = make_T(Mat_MDH)
%MAKE_T 根据改进 DH（MDH）表生成相邻坐标系的符号齐次变换矩阵。
%   输入 Mat_MDH 为 N×4 的符号或数值矩阵，每行依次为：
%       [alpha, a, d, theta]，角度单位为弧度。
%   输出 T 为 4×4×N 的符号数组，T(:,:,i) 对应第 i 行参数。
%   按坐标系从 1 开始编号，T(:,:,i) 为 ^i T_{i+1}，满足：
%       p_i = T(:,:,i) * p_{i+1}（p 为齐次坐标）。
%   使用的 MDH 变换顺序为 Rx(alpha)*Tx(a)*Rz(theta)*Tz(d)。
%
%   示例：
%       T = make_T(Mat_MDH);
%       T_1 = T(:,:,1);  % ^1 T_2
%       T_2 = T(:,:,2);  % ^2 T_3

if ~(isa(Mat_MDH, 'sym') || isnumeric(Mat_MDH)) || ...
        ~ismatrix(Mat_MDH) || size(Mat_MDH, 2) ~= 4 || ...
        size(Mat_MDH, 1) == 0
    error('make_T:InvalidMDH', ...
        'Mat_MDH 必须是非空 N×4 符号或数值矩阵，列顺序为 [alpha, a, d, theta]。');
end

Mat_MDH = sym(Mat_MDH);
N = size(Mat_MDH, 1);
T = sym(zeros(4, 4, N));

for i = 1:N
    alpha = Mat_MDH(i, 1);
    a = Mat_MDH(i, 2);
    d = Mat_MDH(i, 3);
    theta = Mat_MDH(i, 4);

    ca = cos(alpha);
    sa = sin(alpha);
    ct = cos(theta);
    st = sin(theta);

    T(:,:,i) = simplify([ ...
        ct,    -st,     0,    a; ...
        ca*st, ca*ct, -sa, -d*sa; ...
        sa*st, sa*ct,  ca,  d*ca; ...
        0,      0,     0,    1]);
end
end

