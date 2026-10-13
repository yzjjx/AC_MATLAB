function vxf = cross_star(v, f)
%CROSS_STAR 计算空间运动向量与空间力向量的对偶叉乘 v ×* f。
%   输入 v = [omega; u]、f = [力矩; 力]，均为 6×1 符号或数值列向量。
%   输出 vxf 为 6×1 列向量，采用公式2.32：
%       v ×* = [omega_cross, u_cross; zeros(3), omega_cross]
%            = -(v ×).'
%       vxf = (v ×*) * f。
%
%   调用示例：
%       bias_force = cross_star(mat_v{i}, mat_I(:,:,i)*mat_v{i});

if ~(isa(v, 'sym') || isnumeric(v)) || ~isequal(size(v), [6, 1]) || ...
        ~(isa(f, 'sym') || isnumeric(f)) || ~isequal(size(f), [6, 1])
    error('cross_star:InvalidInput', ...
        'v 和 f 必须均为 6×1 的符号或数值列向量。');
end

omega_cross = cross_3x1(v(1:3));
u_cross = cross_3x1(v(4:6));
v_star = [omega_cross, u_cross; zeros(3), omega_cross];
vxf = v_star * f;
end

