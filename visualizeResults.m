function visualizeResults(nodes, elements, U, materialProps)
    numNodes = size(nodes, 1);

    % 提取位移场
    UX = U(1:2:end); % x方向位移
    UY = U(2:2:end); % y方向位移
    dispMag = sqrt(UX.^2 + UY.^2); % 位移幅值

    % 1. 绘制未变形网格和变形网格
    figure;

    % 子图1: 未变形网格
    subplot(1, 2, 1);
    trisurf(elements, nodes(:, 1), nodes(:, 2), zeros(numNodes, 1), 'FaceColor', 'cyan', 'EdgeColor', 'black');
    title('未变形网格');
    xlabel('X方向（单位：米）');
    ylabel('Y方向（单位：米）');
    view(2);
    axis equal;

    % 子图2: 变形网格
    subplot(1, 2, 2);
    deformationScale = 10; % 放大比例
    trisurf(elements, nodes(:, 1) + deformationScale * UX, nodes(:, 2) + deformationScale * UY, ...
            zeros(numNodes, 1), dispMag, 'EdgeColor', 'none');
    colorbar; % 添加色条
    title('变形后的网格（位移幅值）');
    xlabel('X方向（单位：米）');
    ylabel('Y方向（单位：米）');
    view(2);
    axis equal;

    % 2. 绘制应力场
    numElements = size(elements, 1);
    stress = zeros(numElements, 3);
    strain = zeros(numElements, 3);
    D = getMaterialMatrix(materialProps);

    for e = 1:numElements
        elementNodes = elements(e, :);
        coords = nodes(elementNodes, :);
        Ue = [U(2 * elementNodes - 1); U(2 * elementNodes)];
        [B, ~] = calculateBMatrix(coords, [1/3, 1/3]);
        strain(e, :) = B * Ue;
        stress(e, :) = D * strain(e, :)';
    end

    % 绘制应力场 (σx 分量)
    figure;
    trisurf(elements, nodes(:, 1), nodes(:, 2), zeros(numNodes, 1), stress(:, 1), 'EdgeColor', 'none');
    colorbar; % 添加色条
    title('应力场 (\sigma_x)');
    xlabel('X方向（单位：米）');
    ylabel('Y方向（单位：米）');
    view(2);
    axis equal;
end
