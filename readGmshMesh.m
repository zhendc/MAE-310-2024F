function [nodes, elements] = readGmshMesh(mshFile)
    % 读取 Gmsh 生成的 .msh 文件并解析节点和三角形单元信息
    % 输入:
    %   mshFile - Gmsh 生成的 .msh 文件路径
    % 输出:
    %   nodes - Nx3 的矩阵，存储每个节点的坐标 (x, y, z)
    %   elements - Mx3 的矩阵，存储每个三角形单元的节点编号
    % 打开 .msh 文件
    fid = fopen(mshFile, 'r');
    if fid < 0
        error(['无法打开文件: ', mshFile]);
    end
    % 初始化变量
    nodes = [];
    elements = [];
    % 逐行读取文件
    while ~feof(fid)
        % 获取当前行并去除多余空格
        line = strtrim(fgetl(fid));
        % 根据段落标签解析
        switch line
            case '$Nodes' % 解析节点部分
                % 读取节点数量
                numNodes = str2double(fgetl(fid));
                disp(['解析的节点数量为: ', num2str(numNodes)]);
                
                % 初始化节点数组
                nodes = zeros(numNodes, 3);
                for i = 1:numNodes
                    % 读取每个节点的数据
                    data = str2num(fgetl(fid));
                    nodes(i, :) = data(2:4); % 存储 x, y, z 坐标
                end
                fgetl(fid); % 跳过 $EndNodes
            
            case '$Elements' % 解析单元部分
                % 读取单元数量
                numElements = str2double(fgetl(fid));
                disp(['解析的单元数量为: ', num2str(numElements)]);
                
                % 初始化三角形单元数组
                triElements = [];
                
                % 循环读取每个单元
                for i = 1:numElements
                    % 读取一行数据
                    data = str2num(fgetl(fid));
                    
                    % 提取单元类型
                    elementType = data(2); % 第二列是单元类型
                    
                    % 只处理三角形单元（类型 2）
                    if elementType == 2
                        % 三角形单元的最后 3 列是连接的节点编号
                        connectivity = data(end-2:end);
                        triElements = [triElements; connectivity];
                    else
                        % 跳过其他单元类型
                        continue;
                    end
                end
                
                % 存储三角形单元
                elements = triElements;
                fgetl(fid); % 跳过 $EndElements
        end
    end
    
    % 关闭文件
    fclose(fid);
end
% 读取网格文件

