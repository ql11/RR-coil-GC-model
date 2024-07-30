%% 运算环境
    disp('运算环境处理...')
    programdate = datestr(datetime,'yyyymmdd_hhMMss');
    this_path = '/public3/home/scb9962/Upload/RR-coil-GC-model';
    disp(['计算路径：',this_path]);
    addpath(this_path);
    disp('运算环境处理完毕...')


%% 导入数据
disp('正在载入数据...')
    % 计算基础属性
    [Nd,N,Ndp,Nc]= Attitude('N of divisions','SP N','N of DPs','N of coils'); %分解元素数，单饼匝数，双饼个数，线圈个数
    
    % 计算规模提取矩阵选取1号双饼的A单饼和8号双饼的A单饼进行示例计算
    A1 = Nd*N; % 精细模型电流数量也是环向电流数量
    A2 = N; % 粗糙模型电流数量
    
    % 互感矩阵数据
    load([this_path,'/data/mutual_inductance_matrix.mat']);
    % 径向电阻数据
    %load('径向电阻数据')

    % 根据基础属性提取矩阵选取1号双饼的上半单饼和8号双饼的上半单饼进行示例计算
    M_1_1 = M_s_u(1:A1,1:A1);
    M_2_2 = M_s_u(end-N:end,end-N:end);
    M_1_2 = M_s_u(1:A1,end-N:end);
    M_2_1 = M_1_2';

    % 提取径向电阻
    Rr_0 = Attitude('Basic radial resistance');
    R1 = 1/8.*ones(A1,1).*Rr_0; % 径向电阻列向量
    R2 = ones(A2,1).*Rr_0; % 径向电阻列向量

%% 初始化
disp('变量初始化中...')
    % 初始化时间
    t = 0;
    dt = Attitude('Default time step');

    % 初始化电流列向量
    I = zeros(A1+(Nd*(N-1))+A2,1); %电流总数量，表示为α 环向电流-α 径向电流-β 环向电流 

    I(1:A1) = fun_Current_Control(t);% α 环向电流初始值
    I(end-N:end) = fun_Current_Control(t);% β 环向电流初始值


    % 初始化循环数
    loop_number = 100;


%% 数据记录
disp('展开记录文件...')
    % 记录值
    current_record = zeros(size(I,1),loop_number);%电流记录
    time_record = zeros(1,loop_number);%时间轴

%% 循环计算
disp('开始循环运算...')
for loop_count = 1:loop_number

    % 循环时间
        tic
        t = t + dt;
        time_record(loop_count) = t; %记录时间轴
        now_time = datestr(now,'日期yyyy-mm-dd 时间HH:MM:SS');
        disp(['磁体时间：',num2str(t),'秒','本次循环开始时间：【',now_time,' 】']);% 显示当前时间
        prepare_time = toc;tic % 记录准备时间

    % 方程分为9个模块
        % [1,1]
        [equation_1_1,constant_1_1] = level_1_equation(M_1_1,R1);
        % equation_1_1 = equation_1_1 + diag(Circle_Resistance_Matrix_1); 环向电阻-临界电流相关

        % [1,2]

        % [2,1]

        % [2,2]
        equation_2_2 = IC.*M_2_2./dt +...% 一阶微分
        diag(R_2);% 电阻径向电阻
        canstant_2_2 = R_2'*fun_Current_Control(t) +... % 径向电阻乘控制电流
        (IC.*M_2_2*I(N_3 + 1 :N_4))./dt; % 一阶微分
    
        % 组合成为一个方程组
    Equation = [equation_1_1,equation_1_2;...
        equation_2_1,equation_2_2];
        Canstant = [constant_1_1 + canstant_1_2;...
        canstant_2_1 + canstant_2_2];

    % 解方程
    organize_time = toc;tic
    KI = real(Equation\Canstant);

    % 记录  
    current_record(:,loop_count) = KI;%记录电流

    % 显示
    solve_time = toc;tic
    disp(['准备时间',num2str(prepare_time),'秒，构成方程时间',num2str(organize_time),'秒，解方程时间',num2str(solve_time),'秒，本次循环时间',num2str(prepare_time + organize_time + solve_time),'秒']);

end

%% 保存记录
save([this_path,'/record/','record_',programdate,'.mat'],'time_record','current_record');