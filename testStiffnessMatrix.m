function testStiffnessMatrix()
    % 定义节点和单元
    nodes = [0, 0; 1, 0; 0, 1];
    elements = [1, 2, 3];
    materialProps.E = 210e9;
    materialProps.nu = 0.3;
    materialProps.planeStress = true;

    % 测试刚度矩阵
    K_global = assembleStiffnessMatrix(nodes, elements, materialProps);
    assert(isequal(K_global, K_global'), '刚度矩阵不是对称的');

    % 测试特征值正定性
    eigenvalues = eig(K_global);
    assert(all(eigenvalues > 0), '刚度矩阵不是正定的');

    % 测试 Jacobian 和 B 矩阵
    coords = nodes(elements(1, :), :);
    [B, detJ] = calculateBMatrix(coords, [1/3, 1/3]);
    assert(detJ > 0, 'Jacobian 行列式不为正');

    disp('所有测试通过！');
end
