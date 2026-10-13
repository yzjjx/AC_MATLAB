function vxw = cross_6x1(v, w)
%CROSS_6X1 计算两个 6×1 空间运动向量的叉乘。
%   输入 v = [omega; u]、w = [omega_prime; u_prime]，
%   均为 6×1 的符号或数值列向量，前三行为角分量，后三行为线分量。
%   输出 vxw = v × w 为 6×1 列向量，定义为：
%       [omega × omega_prime;
%        omega × u_prime + u × omega_prime]。
%   此函数用于空间运动向量之间的叉乘。

if ~(isa(v, 'sym') || isnumeric(v)) || ~isequal(size(v), [6, 1]) || ...
        ~(isa(w, 'sym') || isnumeric(w)) || ~isequal(size(w), [6, 1])
    error('cross_6x1:InvalidInput', ...
        'v 和 w 必须均为 6×1 的符号或数值列向量。');
end

omega = v(1:3);
u = v(4:6);
omega_prime = w(1:3);
u_prime = w(4:6);

vxw = [cross(omega, omega_prime); ...
       cross(omega, u_prime) + cross(u, omega_prime)];
end

