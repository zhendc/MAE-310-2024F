function [L2_error, H1_error] = calculateErrorNorms(nodes, elements, U, exactSolution, exactGradient, materialProps)
    % 参数:
    %   nodes: 节点坐标 (Nx2)
    %   elements: 单元连接关系 (Mx3)
    %   U: FEM 解的位移场 (2*Nx1)
    %   exactSolution: 精确解函数句柄 @(x, y)
    %   exactGradient: 精确梯度函数句柄 @(x, y) -> [du_dx, du_dy; dv_dx, dv_dy]
    %   materialProps: 材料属性 (用于高斯积分)
    %
    % 输出:
    %   L2_error: L2 范数误差
    %   H1_error: H1 范数误差

    numElements = size(elements, 1);
    L2_error_squared = 0;
    H1_error_squared = 0;
    L2_exact_squared = 0;
    H1_exact_squared = 0;

    [gp, w] = gaussPoints(); % 获取高斯积分点和权重

    for e = 1:numElements
        elementNodes = elements(e, :); % 当前单元的节点编号
        coords = nodes(elementNodes, :); % 当前单元的节点坐标
        Ue = [U(2 * elementNodes - 1); U(2 * elementNodes)]; % 单元位移场

        for i = 1:size(gp, 1)
            % 当前积分点的局部坐标
            xi_eta = gp(i, :);
            weight = w(i);

            % 形状函数计算
            [B, detJ] = calculateBMatrix(coords, xi_eta);

            % 当前积分点的全局坐标
            N = [(1 - sum(xi_eta)), xi_eta(1), xi_eta(2)];
            xy = N * coords; % (x, y)

            % 精确解和梯度
            u_exact = exactSolution(xy(1), xy(2)); % 大小为 1x2
            grad_exact_full = exactGradient(xy(1), xy(2)); % 大小为 2x2
            grad_exact = [grad_exact_full(1, 1); grad_exact_full(2, 2); grad_exact_full(1, 2) + grad_exact_full(2, 1)]; % 转换为应变分量

            % 数值解和梯度
            u_numeric = N * reshape(Ue, 2, [])'; % 修复后的计算
            grad_numeric = B * Ue; % FEM 梯度

            % L2 和 H1 误差累加
            L2_error_squared = L2_error_squared + detJ * weight * sum((u_numeric - u_exact).^2);
            H1_error_squared = H1_error_squared + detJ * weight * (sum((u_numeric - u_exact).^2) + sum((grad_numeric - grad_exact).^2));

            % 精确解的模量
            L2_exact_squared = L2_exact_squared + detJ * weight * sum(u_exact.^2);
            H1_exact_squared = H1_exact_squared + detJ * weight * (sum(u_exact.^2) + sum(grad_exact.^2));
        end
    end

    % 计算相对误差
    L2_error = sqrt(L2_error_squared) / sqrt(L2_exact_squared);
    H1_error = sqrt(H1_error_squared) / sqrt(H1_exact_squared);
end
