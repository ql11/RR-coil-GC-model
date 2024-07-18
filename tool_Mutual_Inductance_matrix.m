% 构建涞水磁体的互感矩阵
% 按照每匝分为8个元素进行互感计算

%% 调用属性参数和初始化
[lx,ly] = Attitude('Length x-axis','Length y-axis'); %x，y轴直线段长度
[r1,dr] = Attitude('Inner fillet radius','Thickness per turn');%圆角内径和每匝厚度

[Nd,N,Ndp,Nc]= Attitude('N of divisions','SP N','N of DPs','N of coils'); %单饼匝数，双饼个数，线圈个数
%N = 10; % 调试用
%[Nd,~,Ndp,Nc]= Attitude('N of divisions','SP N','N of DPs','N of coils'); %单饼匝数，双饼个数，线圈个数

pt = Attitude('Pole pitch');%极距

M_size = N*Nd*(2*Ndp)*Nc;% 矩阵大小，双饼拆为两个单饼
M = zeros(M_size); %初始化互感矩阵,顺序为：元素(1~8)-匝(1~406),矩阵初始位置为位于-x的直线段垂直于x轴
line_array = [1,3,5,7]; %直线段
arc_array = [2,4,6,8]; %圆弧段

% 输出当前日期时间
startDateTime = datetime("now");
fprintf('计算开始 %s\n', startDateTime);

%% 循环计算

for i = 1:M_size % 源元素，矩阵行
%for i = 1:M_size % 源元素，矩阵行
    tic
    %源元素特征
    [Nd_i,N_i,Nsp_i,Ndp_i,Nc_i] = ind2sub([Nd,N,2,Ndp,Nc],i); % 源-方位，匝数，单饼ab面，双饼数，线圈NS极
    r_i = (N_i-1).*dr + r1; %源元素半径,也和直线段位置相关
    pt_i = pt.*(2.*Nc_i-3)/2; %由于线圈NS极导致的位置x坐标变化 N极 -1/2pt S极 1/2pt
    p_i = [lx/2,ly/2].*([-1,1].*(Nd_i==2) + [1,1].*(Nd_i==4) + [1,-1].*(Nd_i==6) + [-1,-1].*(Nd_i == 8)) + [pt_i,0];%圆角圆心位置
    h_i = fun_single_pancake_position(Nsp_i,Ndp_i,Nc_i);
    px_i = fun_Straight_segment_position(Nd_i,Nc_i); % 直线段在坐标轴上的位置
    parfor j = 1:i % 目标元素，矩阵列
        
        %目标元素特征
        [Nd_j,N_j,Nsp_j,Ndp_j,Nc_j] = ind2sub([Nd,N,2,Ndp,Nc],j); % 目标-方位，匝数，单饼ab面，双饼数，线圈NS极
        r_j = (N_j-1).*dr + r1; %目标元素半径，也和直线段位置相关
        pt_j = pt.*(2.*Nc_j-3)/2; %由于线圈NS极导致的位置x坐标变化 N极 -1/2pt S极 1/2pt
        p_j = [lx/2,ly/2].*([-1,1].*(Nd_j==2) + [1,1].*(Nd_j==4) + [1,-1].*(Nd_j==6) + [-1,-1].*(Nd_j == 8)) + [pt_j,0];%圆角圆心位置
        h_j = fun_single_pancake_position(Nsp_j,Ndp_j,Nc_j);
        px_j = fun_Straight_segment_position(Nd_j,Nc_j); % 直线段在坐标轴上的位置
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
                    M(i,j) = 0; % 互感近似为0
                elseif ismember(Nd_i,[1,5])  % 竖直平行
                    l_temp = ly; % 根据方位判断直线段长度
                    
                    d_temp =...
                        abs(r_i - r_j).*(Nd_i == Nd_j) + abs(r_i + r_j).*(Nd_i ~= Nd_j) + ... % 位置差异
                        lx.*(Nd_i ~= Nd_j); % 判断是否位于两侧

                    para_dst = pt.*abs(Nc_i - Nc_j); % 如果不位于同一线圈，则间距为极距
                    M(i,j) = fun_Straight_segment_Mutual_inductance(l_temp,sqrt(d_temp.^2 + (h_i - h_j).^2),para_dst); %平行线互感，叠加高度差
                
                else % 水平平行
                    l_temp = lx; % 根据方位判断直线段长度
                    M(i,j) = fun_Straight_segment_Mutual_inductance(l_temp,sqrt((px_i - px_j).^2 +  (h_i - h_j).^2)); % 平行线互感，间距直接相减，叠加高度差
                end
                
            elseif ismember(Nd_i,arc_array) && ismember(Nd_j,arc_array) % 圆弧段-圆弧段
                afa_i = (Nd_i-2)./2.*pi./2; % 源元素起始角度
                afa_j = (Nd_j-2)./2.*pi./2; % 目标元素起始角度
                p_temp = p_j - p_i;
                M(i,j) = fun_Arc_segment_Mutual_inductance(r_i,r_j,afa_i,afa_i + pi/2,afa_j,afa_j + pi/2,p_temp(1),p_temp(2),h_i - h_j);
            elseif ismember(Nd_i,line_array) &&  ismember(Nd_j,arc_array)% 直线段-圆弧段
                if ismember(Nd_i,[1,5])  % 竖直直线
                    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
                    d_temp = ...
                        r_i - (2*r_i + lx).*... % 根据直线段位置判断位于x、y轴
                        ((abs(Nd_i - Nd_j) - 1) ~= 0) +...% 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_i
                        pt.*(Nc_j - Nc_i); % 竖直直线的极距加到d上
                        
                    M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(ly,r_j,d_temp,h_j-h_i);
                else  % 水平直线                   
                    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
                    d_temp = ...
                        r_i - (2*r_i + ly).*... % 根据直线段位置判断位于x、y轴
                        ((abs(Nd_i - Nd_j) - 1) ~= 0); % 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_i
                    M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(lx,r_j,d_temp,h_j-h_i,(Nc_j - Nc_i));% 竖直直线的极距有额外参数
                end

            elseif ismember(Nd_i,arc_array) &&  ismember(Nd_j,line_array) % 圆弧段-直线段
                if ismember(Nd_j,[1,5])  % 竖直直线
                    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
                    d_temp = ...
                        r_j - (2*r_j + lx).*... % 根据直线段位置判断位于x、y轴
                        ((abs(Nd_j - Nd_i) - 1) ~= 0) +...% 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_j
                        pt.*(Nc_i - Nc_j); % 竖直直线的极距加到d上
                        
                    M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(ly,r_i,d_temp,h_i-h_j);
                else  % 水平直线                   
                    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
                    d_temp = ...
                        r_j - (2*r_j + ly).*... % 根据直线段位置判断位于x、y轴
                        ((abs(Nd_j - Nd_i) - 1) ~= 0); % 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_j
                    M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(lx,r_i,d_temp,h_i-h_j,(Nc_i - Nc_j));% 竖直直线的极距有额外参数
                end

            end
            %M(j,i) = M(i,j);
        end
        
    end
    % 输出进度
    schedule = i/M_size*100;
    fprintf('循环时间%.2fs,循环数%d/%d,进度%.2f%%\n',toc,i,M_size,schedule);
end


%% 保存文件
M_all = tril(M,-1)+M'; % 由下对角线对称镜像形成对称互感矩阵
save([pwd,'\data\mutual_inductance_matrix.mat'],'M_all');
% 输出当前日期时间
endDateTime = datetime('now');
fprintf('计算结束 %s\n', endDateTime);
fprintf('总时长 %s\n', endDateTime-startDateTime);
fprintf('总自感 %.4f\n', sum(M_sp,'all'));
