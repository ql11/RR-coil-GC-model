% 由于`tool_Mutual_Inductance_matrix.m`工具程序生成的`.mat`互感矩阵文件太大（13G）
% 将其进行重构以减小大小，降低后续时间迭代计算的难度
% 重构方式为将3~8双饼的匝内互感划分合并，构成第二等级的电路模型

%% 初始化导入数据
this_path = '/public3/home/scb9962/Upload/RR-coil-GC-model';
disp(['计算路径：',this_path]);
addpath(this_path);

[Nd,N,Ndp,Nc]= Attitude('N of divisions','SP N','N of DPs','N of coils'); %分解元素数，单饼匝数，双饼个数，线圈个数
load([this_path,'/data/mutual_inductance_matrix.mat']);

%% 拆解原互感矩阵
Dn = Nd*N*2*2; % 在此之前的匝每匝分解为8个元素，在此之后的匝每匝不分解
%   1-1 1-2
%   2-1 2-2
M_size = size(M_all);
Dm = (M_size(2) - Dn)/Nd; % 简化后的被简化部分数量

%% 纵向相加 2-1
M_1 = M_all(Dn+1 : end,1:Dn); % 2-1
blockSize = [Nd,size(M_1,2)];
N_1 = blockproc(M_1,blockSize, @(block) sum(block.data,1));

%% 横向相加 1-2
M_2 = M_all(1:Dn,Dn+1 : end); % 1-2
blockSize = [size(M_2,1),Nd];
N_2 = blockproc(M_2,blockSize, @(block) sum(block.data,2));


%% 更简化，将每匝视为一个元素
blockSize = [Nd,Nd];
M_simp = blockproc(M_all,blockSize, @(block) sum(block.data(:)));
%save([this_path,'/data/mutual_inductance_matrix_simple.mat'],'M_simp');

%% 输出
M_mix = [M_all(1:Dn,1:Dn),N_2;N_1,M_simp(Dn/Nd +1 :end,Dn/Nd +1 :end)];
save([this_path,'/data/mutual_inductance_matrix_mix.mat'],'M_mix','-v7.3');

M_mix_single = single(M_mix);
save([this_path,'/data/mutual_inductance_matrix_mix_single_73.mat'],'M_mix_single','-v7.3');
save([this_path,'/data/mutual_inductance_matrix_mix_single.mat'],'M_mix_single');

%% triu
M_s_u = triu(M_mix_single);
save([this_path,'/data/mutual_inductance_matrix_mix_single_triu.mat'],'M_s_u');

%% zip
Msu_str = num2str(M_s_u, '%.6g');
fid = fopen([this_path,'/data/mutual_inductance_matrix_mix_single_str.txt'],'w');
fprintf(fid,'%s\n',Msu_str);
fclose(fid);
gzip([this_path,'/data/mutual_inductance_matrix_mix_single_str.txt']);
