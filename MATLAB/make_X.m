function X = make_X(T)
%MAKE_X 将齐次变换矩阵转换为 Featherstone 空间运动变换矩阵。
%   输入 T 为 4×4 或 4×4×N 的符号或数值数组，每一页为：
%       T(:,:,i) = [R, p; 0, 0, 0, 1] = ^i T_{i+1}。
%   输出 X 为 6×6×N 的符号数组，每一页为：
%       X(:,:,i) = [R.', zeros(3); -R.'*[p]_x, R.']。
%
%   空间运动向量按 [角速度; 线速度] 排列，变换方向为：
%       v_{i+1} = X(:,:,i) * v_i。
%   对应的空间力向量按 [力矩; 力] 排列，反向传递为：
%       f_i = X(:,:,i).' * f_{i+1}。
%   [p]_x 为叉乘矩阵，满足 [p]_x*u = cross(p, u)。
%
%   示例：
%       T = make_T(Mat_MDH);
%       X = make_X(T);
%       X_1 = X(:,:,1);
%       X_2 = X(:,:,2);

if ~(isa(T, 'sym') || isnumeric(T)) || ndims(T) > 3 || ...
        size(T, 1) ~= 4 || size(T, 2) ~= 4 || isempty(T)
    error('make_X:InvalidTransform', ...
        'T 必须是非空 4×4 或 4×4×N 的符号或数值齐次变换矩阵数组。');
end

T = sym(T);
N = size(T, 3);
X = sym(zeros(6, 6, N));
Z = sym(zeros(3));

for i = 1:N
    R = T(1:3, 1:3, i);
    p = T(1:3, 4, i);

    p_cross = [ ...
        0,    -p(3),  p(2); ...
        p(3),  0,    -p(1); ...
       -p(2),  p(1),  0];

    Rt = R.';
    X(:,:,i) = simplify([Rt, Z; -Rt*p_cross, Rt]);
end
end
