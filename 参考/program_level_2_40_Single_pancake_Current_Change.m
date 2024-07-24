% program_level_2_20_Single_pancake_Current_Change.m
% 二级模型的电流变化，用来验证多尺度模型，忽略环向电阻
%% 运算环境
    disp('运算环境处理中...')
    programdate = datestr(datetime,'yyyymmdd_hhMMss');
    this_path = '/public1/home/sc40009/QinLang/TME-model';
    addpath(this_path);

%% 基础参数
    [T,A,r0,Ply,W] = single_pancake_parameter(0); %单饼参数
%% 导入数据
    disp('数据导入中...')
    load('Rough_model_parameter.mat');
    load('level_2_2_Mutual_inductance_Matrix.mat');%level_2_2_Mutual_inductance_Matrix 二级模型之间的互感
    load('level_2_Central_Field_Matrix.mat')%Central_Field_Matrix

    load('level_2_1_Mutual_inductance_Matrix_50.mat');%level_2_1_Mutual_inductance_Matrix 二级模型与一级模型之间的互感
    load('level_2_2_Vertical_Field.mat');%level_2_2_Vertical_Field 二级模型造成的垂直磁场
    load('level_2_thermal_conductivity_matrix.mat');%level_2_thermal_conductivity_matrix 二级模型热传导矩阵
    disp('数据导入完毕！')
%% 初始赋值
    disp('数据初始化中...')
    %时间初始值
    t = 0;

    %循环值
    loop_number = 1440; %循环次数，可以估计运算时间

    %α面下面的β面个数
    Beta_up  = 0;
    Beta_down = 40; 

    N = (Beta_up + Beta_down)*T; 
    [Mutual_Inductance_Matrix_2_1,Mutual_Inductance_Matrix_2_2] = fun_Mutual_Inductance_beta(level_2_1_Mutual_inductance_Matrix,level_2_2_Mutual_inductance_Matrix,Rough_model_Mutual_inductance_Matrix,Beta_up,Beta_down);% 二级对一级电感 二级之间的电感
    % Mutual_Inductance_Matrix_2_1 二级对一级电感
    % Mutual_Inductance_Matrix_2_2 二级之间的电感
    Circle_Resistance_Matrix_2 = zeros(N,1);% β的初始环向电阻

    % 根据一级模型的位置调整二级模型对一级模型的垂直场和互感影响
    Vertical_Field_Matrix_2_1 = fun_VMF_beta(level_2_2_Vertical_Field,Beta_up,Beta_down - 1);% 二级垂直场计算矩阵

    %电流初值
    I = zeros(N,1);% 电流列向量大小

    % 记录值
    current_record = zeros(N,loop_number);%电流记录
    time_record = zeros(1,loop_number);%时间轴
    voltage_record = zeros(1,loop_number);
    field_record = zeros(1,loop_number);
    temperature_record = zeros(N,loop_number);%温度记录
    disp('数据初始化完毕!开始循环运算...')

    % 电压显示
    %fig_V = figure('Name','电压');
%% 循环计算
    for loop_count = 1:loop_number
        tic
        dt = 0.25;%fun_Time_Interval_Control('time',t); %此次循环对应的时间间隔
        t = t + dt;
        time_record(loop_count) = t; %记录时间轴
        now_time = datestr(now,'日期yyyy-mm-dd 时间HH:MM:SS');
        disp(['磁体时间：',num2str(t),'秒','本次循环开始时间：【',now_time,' 】']);
        prepare_time = toc;tic
        
        % 系数包括：互感和自感/dt+径向电阻
        % 常数包括：互感和自感乘上个时间点的电流/dt+控制电流*径向电阻
        % 常数
        Canstant =  Mutual_Inductance_Matrix_2_2*I./dt+...
                    repmat(fun_Current_Control(t).*Rough_model_Radial_resistance,Beta_up + Beta_down,1);
        % 互感和自感/dt
        Equation = Mutual_Inductance_Matrix_2_2./dt;
        % 加上径向电阻
        Equation = Equation + diag(repmat(Rough_model_Radial_resistance,Beta_up + Beta_down,1));

        % 环向电阻需要计算温度和垂直场
        for k = 1:(Beta_up+Beta_down)
            Temperature_temp = fun_level_2_Tamperature_Change(I((k-1)*T+1:k*T),level_2_thermal_conductivity_matrix,Rough_model_Radial_resistance,Circle_Resistance_Matrix_2((k-1)*T+1:k*T),fun_Current_Control(t));
            % 垂直磁场
            VMF_temp = Vertical_Field_Matrix_2_1(:,:,k)*I;% 垂直场，注意不包括单饼的垂直自场
            % 迭代环向电阻 %环形电阻应当是一级模型的形式参与计算 迭代用的环形电阻应当是二级模型的形式
            Circle_Resistance_Matrix_2((k-1)*T+1:k*T) = fun_Circle_Resistance_Matrix(I((k-1)*T + 1 : k*T),VMF_temp,Temperature_temp);
            Temperature((k-1)*T+1:k*T) = Temperature_temp;
        end
        % 在方程里面加上径向电阻
        Equation = Equation + diag(Circle_Resistance_Matrix_2);




        organize_time = toc;tic
        I = real(Equation\Canstant);
        current_record(:,loop_count) = I;%记录电流
        now_voltage = (fun_Current_Control(t) - I)'*repmat(Rough_model_Radial_resistance,(Beta_up + Beta_down),1);
        voltage_record(:,loop_count) = now_voltage;%记录总电压
        now_field = fun_Field_Multi_pancake(I,level_2_Central_Field_Matrix);
        field_record(:,loop_count) = now_field;%记录中心磁场
        temperature_record(:,loop_count) = Temperature;%记录温度

        %if rem(loop_count,10) == 0
            %figure(fig_V);
            %scatter(time_record,voltage_record,'.') 
        %end

        solve_time = toc;tic
        disp(['准备时间',num2str(prepare_time),'秒，构成方程时间',num2str(organize_time),'秒，解方程时间',num2str(solve_time),'秒，本次循环时间',num2str(prepare_time + organize_time + solve_time),'秒']);
    
        disp(['V=',num2str(now_voltage),' | B=',num2str(now_field)])
    end

    save([this_path,'/record/','record_',programdate,'.mat'],'time_record','current_record','voltage_record','field_record','temperature_record');
    
