function K_global = assembleStiffnessMatrix(nodes, elements, materialProps)
    % 参数:
    % nodes: Nx2矩阵，每行是节点的 (x, y) 坐标
    % elements: Mx3矩阵，每行是一个三角形单元的节点编号
    % materialProps: 材料属性结构体 (E, nu, planeStress)
    
    numNodes = size(nodes, 1);  % 节点数量
    numElements = size(elements, 1);  % 单元数量
    K_global = zeros(2 * numNodes, 2 * numNodes);  % 初始化全局刚度矩阵

    D = getMaterialMatrix(materialProps);  % 材料矩阵
    [gp, w] = gaussPoints();  % 高斯积分点和权重

    for e = 1:numElements  % 遍历每个单元
        elementNodes = elements(e, :);  % 当前单元的节点编号
        coords = nodes(elementNodes, :);  % 当前单元的节点坐标
        
        % 初始化局部刚度矩阵
        Ke_local = zeros(6, 6);

        for i = 1:size(gp, 1)  % 遍历每个高斯积分点
            [B, detJ] = calculateBMatrix(coords, gp(i, :));  % 计算 B 矩阵和 Jacobian 行列式
            
            % 调试信息: 打印 Jacobian 行列式
            disp(['单元 ', num2str(e), ', 积分点 ', num2str(i), ': Jacobian 行列式 = ', num2str(detJ)]);
            
            % 验证 Jacobian 行列式是否为正
            if detJ <= 0
                error(['单元 ', num2str(e), ' 的 Jacobian 行列式为负或零，请检查网格和节点顺序！']);
            end
            
            % 累加局部刚度矩阵
            Ke_local = Ke_local + (B' * D * B) * detJ * w(i);
        end

        % 调试信息: 打印局部刚度矩阵
        disp(['单元 ', num2str(e), ' 的局部刚度矩阵:']);
        disp(Ke_local);

        % 将局部刚度矩阵组装到全局刚度矩阵
        dof = [elementNodes * 2 - 1; elementNodes * 2];  % 自由度索引
        dof = dof(:);  % 转换为列向量
        K_global(dof, dof) = K_global(dof, dof) + Ke_local;
    end
end

