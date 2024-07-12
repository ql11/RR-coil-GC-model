% 构建涞水磁体单饼的互感矩阵
% 按照每匝分为8个元素进行互感计算

%% 调用属性参数和初始化
[lx,ly] = Attitude('Length x-axis','Length y-axis'); %x，y轴直线段长度
[r1,dr] = Attitude('Inner fillet radius','Thickness per turn');%圆角内径和每匝厚度
[Nd,N]= Attitude('N of divisions','SP N'); %单饼匝数，双饼个数，线圈个数

M_size = N*Nd;% 矩阵大小，双饼拆为两个单饼
M = zeros(M_size); %初始化互感矩阵,顺序为：元素(1~8)-匝(1~406),矩阵初始位置为位于-x的直线段垂直于x轴

%% 循环计算

for i = 1:M_size % 源元素，矩阵行
    
    %源元素特征
    [Nd_i,N_i] = ind2sub([Nd,N],i); % 源方位，源匝数
    r_i = (N_i-1).*dr + r1; %源元素半径
    
    for j = 1:i % 目标元素，矩阵列
        
        %目标元素特征
        [Nd_j,N_j] = ind2sub([Nd,N],j); % 目标方位，目标匝数
        r_j = (N_j-1).*dr + r1; %目标元素半径
        
        if i == j
            if ismember(Nd_i,[1,5]) %直线段 y方向
                M(i,j) = fun_Straight_segment_Mutual_inductance(ly);
            elseif ismember(Nd_i,[3,7]) %直线段 x方向
                M(i,j) = fun_Straight_segment_Mutual_inductance(lx);
            else  % 圆弧段
                M(i,j) = fun_Arc_segment_Mutual_inductance(r_i,r_j);
            end
        else
            if  ismember(Nd_i,[1,3,5,7]) && ismember(Nd_j,[1,3,5,7]) %均为直线
                if mod(Nd_i - Nd_j, 4) ~= 0 % 垂直
                    M(i,j) = 0; % 互感为0
                else % 平行
                    l_temp = lx.*(ismember(Nd_i,[3,7])) + ly.*(ismember(Nd_i,[1,5])); % 根据方位判断直线段长度
                    d_temp = abs(r_i-r_j) + ... % 位置差异
                    (lx.*(ismember(Nd_i,[1,5])) + ly.*(ismember(Nd_i,[3,7])))... % xy方位判断
                    .*abs(mod(Nd_i - Nd_j,8)/4); % 判断是否位于两侧
                    M(i,j) = fun_Straight_segment_Mutual_inductance(l_temp,d_temp); %平行线互感
                end
%             elseif condition2 % 平行直线
%                 body
%             elseif condition3 % 圆弧段-圆弧段
%                 body
%             elseif condition4 % 直线段-圆弧段
%                 body
%             elseif condition5 % 圆弧段-直线段
%                 body
            end
            M(j,i) = M(i,j);
        end
        
    end
end