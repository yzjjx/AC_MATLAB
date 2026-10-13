function c_cross = cross_3x1(c)
%CROSS_3X1 将 3×1 向量转换为 3×3 反对称叉乘矩阵。
%   输入 c = [cx; cy; cz]，为符号或数值列向量。
%   输出 c_cross = [0, -cz, cy; cz, 0, -cx; -cy, cx, 0]，
%   满足 c_cross*u = cross(c, u)。
%
%   调用示例：
%       c_cross = cross_3x1(mat_c(:,i));

if ~(isa(c, 'sym') || isnumeric(c)) || ~isequal(size(c), [3, 1])
    error('cross_3x1:InvalidInput', ...
        'c 必须是 3×1 的符号或数值列向量，例如 mat_c(:,i)。');
end

c_cross = [ ...
     0,    -c(3),  c(2); ...
     c(3),  0,    -c(1); ...
    -c(2),  c(1),  0];
end

