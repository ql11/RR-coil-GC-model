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
    [N,Ndp,Nc]= Attitude('SP N','N of DPs','N of coils'); %分解元素数，单饼匝数，双饼个数，线圈个数
    
    % 计算规模提取矩阵选取1号双饼的A单饼和8号双饼的A单饼进行示例计算
    A = N*Ndp*2*Nc; % 精细模型电流数量也是环向电流数量
    
    % 互感矩阵数据
    load([this_path,'/data/mutual_inductance_matrix_simple.mat']);
    M = M_simp;

    % 提取径向电阻
    Rr_0 = Attitude('Basic radial resistance');
    R = ones(A,1).*Rr_0; % 径向电阻列向量

%% 初始化
disp('变量初始化中...')
    % 初始化时间
    t = 0;
    dt = Attitude('Default time step');

    % 初始化电流列向量
    I = oness(A,1).*Operation_Current(t); %环向电流总数量

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

    % 电流接近临界电流时的环向电阻，未定义，用0代替
    Rc = zeros(size(R));

    % 方程系数
    Equation = M./dt + diag(Rc) + diag(R); % 互感电压加失超电压；径向电阻项由Canstant移项而来，无独立物理意义
    
    % 方程常数
    Canstant = M*I./dt + Operation_Current(t).*R;

    % 解方程
    organize_time = toc;tic
    KI = real(Equation\Canstant);

    % 记录  
    current_record(:,loop_count) = KI;%记录电流
    I = KI;

    % 显示
    solve_time = toc;tic
    disp(['准备时间',num2str(prepare_time),'秒，构成方程时间',num2str(organize_time),'秒，解方程时间',num2str(solve_time),'秒，本次循环时间',num2str(prepare_time + organize_time + solve_time),'秒']);

end

%% 保存记录
save([this_path,'/record/','record_',programdate,'.mat'],'time_record','current_record');