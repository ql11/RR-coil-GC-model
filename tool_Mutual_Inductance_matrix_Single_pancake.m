% 构建涞水磁体单饼的互感矩阵
% 按照每匝分为8个元素进行互感计算

%% 调用属性参数和初始化
[lx,ly] = Attitude('Length x-axis','Length y-axis'); %x，y轴直线段长度
[r1,dr] = Attitude('Inner fillet radius','Thickness per turn');%圆角内径和每匝厚度
[Nd,N]= Attitude('N of divisions','SP N'); %单饼匝数，双饼个数，线圈个数

N = 10; % 调试用

M_size = N*Nd;% 矩阵大小，双饼拆为两个单饼
M = zeros(M_size); %初始化互感矩阵,顺序为：元素(1~8)-匝(1~406),矩阵初始位置为位于-x的直线段垂直于x轴
line_array = [1,3,5,7]; %直线段
arc_array = [2,4,6,8]; %圆弧段


%% 循环计算

for i = 1:M_size % 源元素，矩阵行
    
    %源元素特征
    [Nd_i,N_i] = ind2sub([Nd,N],i); % 源方位，源匝数
    r_i = (N_i-1).*dr + r1; %源元素半径,也和直线段位置相关
    p_i = [lx/2,ly/2].*([-1,1].*(Nd_i==2) + [1,1].*(Nd_i==4) + [1,-1].*(Nd_i==6) + [-1,-1].*(Nd_i == 8));%圆角圆心位置
    
    for j = 1:i % 目标元素，矩阵列
        
        %目标元素特征
        [Nd_j,N_j] = ind2sub([Nd,N],j); % 目标方位，目标匝数
        r_j = (N_j-1).*dr + r1; %目标元素半径，也和直线段位置相关
        p_j = [lx/2,ly/2].*([-1,1].*(Nd_j==2) + [1,1].*(Nd_j==4) + [1,-1].*(Nd_j==6) + [-1,-1].*(Nd_j == 8));%圆角圆心位置


        if i == j % 对角线计算自感
            if ismember(Nd_i,[1,5]) %直线段 y方向
                M(i,j) = fun_Straight_segment_Mutual_inductance(ly);
            elseif ismember(Nd_i,[3,7]) %直线段 x方向
                M(i,j) = fun_Straight_segment_Mutual_inductance(lx);
            else  % 圆弧段
                M(i,j) = fun_Arc_segment_Mutual_inductance(r_i,r_j);
            end
        else % 非对角线计算互感
            if  ismember(Nd_i,line_array) && ismember(Nd_j,line_array) %均为直线
                if mod(Nd_i - Nd_j, 4) ~= 0 % 垂直直线
                    M(i,j) = 0; % 互感为0
                else % 非垂直直线即为平行直线
                    l_temp = lx.*(ismember(Nd_i,[3,7])) + ly.*(ismember(Nd_i,[1,5])); % 根据方位判断直线段长度
                    
                    d_temp =...
                        abs(r_i-r_j).*(Nd_i == Nd_j) + abs(r_i-r_j).*(Nd_i ~= Nd_j) + ... % 位置差异
                        (lx.*(ismember(Nd_i,[1,5])) + ly.*(ismember(Nd_i,[3,7])))... % xy方位判断
                        .*(Nd_i ~= Nd_j); % 判断是否位于两侧
                    
                    M(i,j) = fun_Straight_segment_Mutual_inductance(l_temp,d_temp); %平行线互感

                    clear l_temp d_temp
                end
            elseif ismember(Nd_i,arc_array) && ismember(Nd_j,arc_array) % 圆弧段-圆弧段
                afa_i = (Nd_i-2)./2.*pi./2; % 源元素起始角度
                afa_j = (Nd_j-2)./2.*pi./2; % 目标元素起始角度
                p_temp = p_j - p_i;
                M(i,j) = fun_Arc_segment_Mutual_inductance(r_i,r_j,afa_i,afa_i + pi/2,afa_j,afa_j + pi/2,p_temp(1),p_temp(2),0);
            elseif ismember(Nd_i,line_array) &&  ismember(Nd_j,arc_array)% 直线段-圆弧段

                l_temp = ly.*ismember(Nd_i,[1,5]) + lx.*ismember(Nd_i,[3,7]); % 直线段长度位于y轴线上采用lx,位于x轴线上采用ly
                    
                    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
                d_temp = ...
                    r_i - (2*r_i + lx.*ismember(Nd_i,[1,5]) + ly.*ismember(Nd_i,[3,7])).*... % 根据直线段位置判断位于x、y轴
                    ((abs(Nd_i - Nd_j) - 1) ~= 0); % 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_i
                    
                M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(l_temp,r_j,d_temp,0);
 
            elseif ismember(Nd_i,arc_array) &&  ismember(Nd_j,line_array) % 圆弧段-直线段
                
                l_temp = ly.*ismember(Nd_j,[1,5]) + lx.*ismember(Nd_j,[3,7]); % 直线段长度位于y轴线上采用lx,位于x轴线上采用ly
                d_temp = ...
                    r_j - (2*r_j + lx.*ismember(Nd_j,[1,5]) + ly.*ismember(Nd_j,[3,7])).*... % 根据直线段位置判断位于x、y轴
                    ((abs(Nd_i - Nd_j) - 1) ~= 0); % 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_j
                
                M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(l_temp,r_i,d_temp,0);
            end
            M(j,i) = M(i,j);
        end
        
    end
end