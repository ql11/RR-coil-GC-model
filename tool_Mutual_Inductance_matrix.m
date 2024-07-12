% 构建涞水磁体的互感矩阵
% 按照每匝分为8个元素进行互感计算

%% 调用属性参数和初始化
[lx,ly] = Attitude('Length x-axis','Length y-axis'); %x，y轴直线段长度
[r1,dr] = Attitude('Inner fillet radius','Thickness per turn');%圆角内径和每匝厚度
[Nd,N,Ndp,Nc]= Attitude('N of divisions','SP N','N of DPs','N of coils'); %单饼匝数，双饼个数，线圈个数
pt = Attitude('Pole pitch');%极距

%todo 其他参数未调用

M_size = N*2*Ndp*Nc*Nd;% 矩阵大小，双饼拆为两个单饼

M = zeros(M_size); %初始化互感矩阵,顺序为：元素(1~8)-匝(1~406)-单饼(1~2)-双饼(1~4)-线圈(1~2)

%% 循环计算

for i = 1:M_size % 原元素，矩阵行
    %互感计算
    [Nd_i,N_i,Nsp_i,Ndp_i,Nc_i] = ind2sub([Nd,N,2,Ndp,Nc],i);
    for j = 1:i % 新元素，矩阵列
        %todo
        [Nd_j,N_j,Nsp_j,Ndp_j,Nc_j] = ind2sub([Nd,N,2,Ndp,Nc],j);
        
        if i == j
            M(i,j) = 0; % 自感
        elseif condition1 %垂直直线
            M(i,j) = 0;
        elseif condition2 % 平行直线
            body
        elseif condition3 % 圆弧段-圆弧段
            body
        elseif condition4 % 直线段-圆弧段
            body
        elseif condition5 % 圆弧段-直线段
            body
        end
        M(j,i) = M(i,j);
    end
end

M = fun_Arc_segment_Mutual_inductance(R1,R2,afa1,afa2,bta1,bta2,dx,dy,dz)