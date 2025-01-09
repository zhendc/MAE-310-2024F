function [K_mod, F_mod] = applyBoundaryConditions(K_global, F_global, boundaryConditions, numNodes)
    % 参数:
    %   K_global - 全局刚度矩阵
    %   F_global - 全局力向量
    %   boundaryConditions - 包含 Dirichlet 和 Neumann 边界条件的结构体
    %   numNodes - 总节点数

    totalDofs = 2 * numNodes;  % 每个节点有两个自由度
    fixedDofs = [];  % 存储被固定的自由度索引

    % 1. 应用 Dirichlet 边界条件
    if isfield(boundaryConditions, 'Dirichlet')
        for i = 1:length(boundaryConditions.Dirichlet)
            node = boundaryConditions.Dirichlet(i).node;
            dof = boundaryConditions.Dirichlet(i).dof;
            fixedDofs = [fixedDofs; 2 * node - (2 - dof)];
        end
    end

    % 计算自由自由度索引
    freeDofs = setdiff(1:totalDofs, fixedDofs);

    % 2. 应用 Neumann 边界条件
    if isfield(boundaryConditions, 'Neumann')
        for i = 1:length(boundaryConditions.Neumann)
            node = boundaryConditions.Neumann(i).node;
            dof = boundaryConditions.Neumann(i).dof;
            value = boundaryConditions.Neumann(i).value;
            F_global(2 * node - (2 - dof)) = F_global(2 * node - (2 - dof)) + value;
        end
    end

    % 3. 修改刚度矩阵和力向量
    K_mod = K_global(freeDofs, freeDofs);
    F_mod = F_global(freeDofs);
end
