% 测试脚本
disp('测试开始...');

% 1. 定义网格数据
disp('加载测试网格...');
nodes = [0, 0; 1, 0; 1, 1; 0, 1]; % 节点坐标
triangleElements = [1, 2, 3; 2, 4, 3]; % 三角形单元
pointElements = [1; 4]; % 点单元（节点编号）
lineElements = [1, 2; 2, 3; 3, 4; 4, 1]; % 线单元（边界）

disp('节点坐标:');
disp(nodes);
disp('三角形单元:');
disp(triangleElements);
disp('点单元:');
disp(pointElements);
disp('线单元:');
disp(lineElements);

% 2. 材料属性
disp('定义材料属性...');
materialProps.E = 210e9; % 弹性模量 (Pa)
materialProps.nu = 0.3; % 泊松比
materialProps.planeStress = true; % 平面应力分析
disp(materialProps);

% 3. 计算全局刚度矩阵
disp('计算全局刚度矩阵...');
K_global = assembleStiffnessMatrix(nodes, triangleElements, materialProps, lineElements, pointElements);
disp('全局刚度矩阵 K_global (前 5x5):');
disp(K_global(1:min(5, size(K_global, 1)), 1:min(5, size(K_global, 2))));

% 4. 初始化全局力向量
disp('初始化全局力向量...');
F_global = zeros(size(K_global, 1), 1);
disp('全局力向量 F_global (前 5):');
disp(F_global(1:min(5, length(F_global))));

% 5. 定义边界条件
disp('定义边界条件...');
boundaryConditions.Dirichlet = [...
    struct('node', 1, 'dof', 1, 'value', 0);  % 固定点单元 1 的 x 自由度
    struct('node', 1, 'dof', 2, 'value', 0);  % 固定点单元 1 的 y 自由度
    struct('node', 4, 'dof', 2, 'value', 0);  % 固定点单元 4 的 y 自由度
];
boundaryConditions.Neumann = [...
    struct('node', 2, 'dof', 2, 'value', -1000); % 在节点 2 上施加 -1000 N 的 y 方向外力
];
disp('Dirichlet 边界条件:');
disp(boundaryConditions.Dirichlet);
disp('Neumann 边界条件:');
disp(boundaryConditions.Neumann);

% 6. 应用边界条件
disp('应用边界条件...');
[K_mod, F_mod] = applyBoundaryConditions(K_global, F_global, boundaryConditions, size(nodes, 1));
disp('修改后的刚度矩阵 K_mod (前 5x5):');
disp(K_mod(1:min(5, size(K_mod, 1)), 1:min(5, size(K_mod, 2))));
disp('修改后的力向量 F_mod (前 5):');
disp(F_mod(1:min(5, length(F_mod))));

% 7. 求解位移场
disp('求解位移场...');
U = K_mod \ F_mod;
disp('位移场 U (前 5):');
disp(U(1:min(5, length(U))));

disp('测试结束。');
