function boundaryConditions = specifyBoundaryConditions()
    % 用户指定边界条件
    fprintf('请指定 Dirichlet 边界条件:\n');
    fprintf('输入格式: [节点编号, 自由度(1=x, 2=y), 位移值]\n');
    fprintf('示例: [1, 1, 0; 1, 2, 0; 3, 2, 0]\n');
    dirichletInput = input('输入 Dirichlet 边界条件 (用矩阵格式): ');
    boundaryConditions.Dirichlet = struct([]);
    for i = 1:size(dirichletInput, 1)
        boundaryConditions.Dirichlet(i).node = dirichletInput(i, 1);
        boundaryConditions.Dirichlet(i).dof = dirichletInput(i, 2);
        boundaryConditions.Dirichlet(i).value = dirichletInput(i, 3);
    end

    fprintf('请指定 Neumann 边界条件:\n');
    fprintf('输入格式: [节点编号, 自由度(1=x, 2=y), 力值]\n');
    fprintf('示例: [2, 1, 1000]\n');
    neumannInput = input('输入 Neumann 边界条件 (用矩阵格式): ');
    boundaryConditions.Neumann = struct([]);
    for i = 1:size(neumannInput, 1)
        boundaryConditions.Neumann(i).node = neumannInput(i, 1);
        boundaryConditions.Neumann(i).dof = neumannInput(i, 2);
        boundaryConditions.Neumann(i).value = neumannInput(i, 3);
    end
end
