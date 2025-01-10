function K_global = assembleStiffnessMatrix(nodes, elements, materialProps)
    % 参数:
    %   nodes: Nx2矩阵，每行是节点的 (x, y) 坐标
    %   elements: Mx3矩阵，每行是一个三角形单元的节点编号
    %   materialProps: 材料属性结构体 (E, nu, planeStress)
    %
    % 输出:
    %   K_global: 全局刚度矩阵

    % 节点数量和单元数量
    numNodes = size(nodes, 1);
    numElements = size(elements, 1);

    % 初始化全局刚度矩阵
    K_global = zeros(2 * numNodes, 2 * numNodes);

    % 获取材料矩阵
    D = getMaterialMatrix(materialProps);

    % 获取高斯积分点和权重
    [gp, w] = gaussPoints();

    % 遍历每个单元
    for e = 1:numElements
        % 获取当前单元的节点编号
        elementNodes = elements(e, :);

        % 获取当前单元的节点坐标
        coords = nodes(elementNodes, :);

        % 初始化当前单元的局部刚度矩阵
        Ke_local = zeros(6, 6);

        % 遍历每个高斯积分点
        for i = 1:size(gp, 1)
            % 计算形状函数的 B 矩阵和 Jacobian 行列式
            [B, detJ] = calculateBMatrix(coords, gp(i, :));

            % 检查 Jacobian 行列式是否有效
            if detJ <= 0
                error(['单元 ', num2str(e), ' 的 Jacobian 行列式为负或零，请检查网格和节点顺序！']);
            end

            % 累加到局部刚度矩阵
            Ke_local = Ke_local + (B' * D * B) * detJ * w(i);
        end

        % 将局部刚度矩阵组装到全局刚度矩阵
        dof = [elementNodes * 2 - 1; elementNodes * 2];  % 每个节点的自由度索引
        dof = dof(:);  % 转为列向量
        K_global(dof, dof) = K_global(dof, dof) + Ke_local;
    end
end

