% program_4_level_magnet_Current_Change.m

%% 运算环境
disp('运算环境处理...')
    programdate = datestr(datetime,'yyyymmdd_hhMMss');
    this_path =  '/public1/home/sc40009/QinLang/TME-model'; % A区
    %this_path = '/public1/home/sc40009/QinLang/From Lptp_2/TME-model'; % D区
    addpath(this_path);
    disp('运算环境处理完毕...')
%% 导入数据
disp('正在载入数据...')

    IC = 1./1.09;
    % Lv.1模型
    load('Single_pancake_Mutual_inductance_Matrix.mat','-mat'); %一级模型每个元素的互感矩阵，包括自感
    load('Radial_resistance_Matrix.mat','-mat'); %每个元素的径向电阻
    load('self_Vertical_Magnet_field_Matrix.mat');% 自场垂直分量计算矩阵
    load('thermal_conductivity_matrix.mat');% 热传导矩阵

    % lv.2模型
    load('Rough_model_parameter.mat');
    %couple_inductance_Matrix % 粗糙模型与精细模型之间的互感 [T*A,T]
    %Rough_model_Mutual_inductance_Matrix % 粗糙模型互感与自感
    %Rough_model_Radial_resistance % 粗糙模型的径向电阻
    load('level_2_2_Mutual_inductance_Matrix.mat');%level_2_2_Mutual_inductance_Matrix 二级模型之间的互感
    load('level_2_2_Vertical_Field_Titan.mat');%level_2_2_Vertical_Field 二级模型造成的垂直磁场
    load('level_2_thermal_conductivity_matrix.mat');%level_2_thermal_conductivity_matrix 二级模型热传导矩阵

    load('Double_pancake_Central_Magnet_field_Matrix.mat');%中心磁场计算，既有精细模型也有粗糙模型
    load('level_2_1_Mutual_inductance_Matrix_50.mat');%level_2_1_Mutual_inductance_Matrix 二级模型与一级模型之间的互感
    load('level_2_Central_Field_Matrix.mat')%Central_Field_Matrix

    % lv.3模型
    load('level_3_1_Mutual_inductance_Matrix.mat')
    load('level_3_2_Mutual_inductance_Matrix.mat')
    load('level_3_2_Vertical_field.mat')
    load('level_3_Central_magnetic_field_Matrix.mat') % level_3_Central_magnetic_field_inner (1,8)内层双绕单饼对中心磁场的贡献；level_3_Central_magnetic_field_outter (1,72) 外层线圈对中心磁场的贡献
    load('level_3_Mutual_inductance_Matrix.mat') %level_3_Mutual_inductance_Matrix (80,80)lv.3互感矩阵
    load('level_3_Radial_Resistance_Matrix.mat')
    % level_3_Radial_Resistance_Matrix_inner (1,8)电阻；level_3_Radial_Resistance_Matrix_outter (1,72)电阻
    load('level_2_3_Vertical_field_fix.mat')
    load('level_3_3_Vertical_field.mat')

    % lv.0模型
    load('lv0_Mutual_inductance_Matrix.mat'); % lv0模型互感 lv0_Mutual_inductance_Matrix
    % lv0_Mutual_inductance_Matrix(T*A - A +A*40 + 39*T + 1:T*A - A +A*40 + 39*T + 4,:) = lv0_Mutual_inductance_Matrix(T*A - A +A*40 + 39*T + 1:T*A - A +A*40 + 39*T + 4,:)./10;
    disp('数据导入完毕！')

%% 初始化数值
disp('变量初始化中...')
    [T,A,r0,Ply,W] = single_pancake_parameter(0); %单饼参数
    %时间初始值
    t = 0;
    %循环值
    loop_number = 2000; %循环次数，可以估计运算时间
    %磁体状态
    n_1 = 40;% 内层单绕单饼
    n_2 = 8;% 内层双绕单饼
    n_3 = 72;% 外层单饼
    %α面下面的β面个数
    Beta_up  = 0;
    Beta_down = n_1 - Beta_up - 1; % 39

    % lv0细分数量
    [lv0_p,~,~] = lv0_position();
    lv0_n = 40;

    %常用数值
    N_1 = T*A; % 对应α面 lv.1 模型的环向电流数量
    N_0 = T*A - A + A*lv0_n;%对应α面 lv.0 模型的环向电流
    N_2 = (2*T-1)*A;% 对应 α面 精细模型的 全部电流
    N_3 = N_2 + T*Beta_down;% lv.1 + lv.2电流个数
    N_4 = N_3 + n_2 + n_3;% lv.1 + lv.2 + lv.3电流个数

    P_0 = A*(lv0_p - 1); % lv0 部分开始位置-1
    P_1 = P_0 + A*lv0_n; %lv0部分结束的位置
    P_2 = T*A - A + A*lv0_n; % 全部电流数量  

    %电流初值
    I = zeros(N_4,1);% 电流列向量大小
    I(1:N_1) = fun_Current_Control(t);% α 环向电流初始值
    I(N_2 + 1:N_4) = fun_Current_Control(t);% β 环向电流初始值

    KI = fun_I_format_change_01(I,zeros(N_4 - A + A*lv0_n,1));
    Temperature = ones(N_1,1).*External_temperature(0);%温度初始值，只对应α面的环向电流

    Circle_Resistance_Matrix_1 = fun_Circle_Resistance_Matrix(I(1:N_2),0,Temperature);% α初始环向电阻，认为垂直磁场为0和温度为0外界温度，在第一次循环之后就能够有正确的温度和磁场了TT
    Circle_Resistance_Matrix_2 = repmat(fun_Circle_Resistance_Matrix(I(N_2 + 1:N_2 + T),0,ones(T,1).*External_temperature(0)),Beta_up+Beta_down,1);% β的初始环向电阻

    %% 数据预处理
disp('数据预处理中...')
    % 根据一级模型的位置调整二级模型对一级模型的垂直场和互感影响
    Vertical_Field_Matrix_2_1 = fun_VMF_beta(level_2_2_Vertical_Field,Beta_up,Beta_down);% 二级垂直场计算矩阵

    [Mutual_Inductance_Matrix_2_1,Mutual_Inductance_Matrix_2_2] = fun_Mutual_Inductance_beta(IC.*level_2_1_Mutual_inductance_Matrix,IC.*level_2_2_Mutual_inductance_Matrix,IC.*Rough_model_Mutual_inductance_Matrix,Beta_up,Beta_down);% 二级对一级电感 二级之间的电感
    % Mutual_Inductance_Matrix_2_1 二级对一级电感
    % Mutual_Inductance_Matrix_2_2 二级之间的电感
    M_2_1 = Mutual_Inductance_Matrix_2_1';M_2_0 = [M_2_1(:,1:P_0),lv0_Mutual_inductance_Matrix(T*A+1:T*A+T*Beta_down,:),M_2_1(:,lv0_n*A + A + 1 :end)];M_0_2 = M_2_0';

    level_3_Radial_Resistance_Matrix = [level_3_Radial_Resistance_Matrix_inner,level_3_Radial_Resistance_Matrix_outter]; % 三级模型径向电阻 (1,80)
    
    level_3_Central_magnetic_field_Matrix = [level_3_Central_magnetic_field_inner,level_3_Central_magnetic_field_outter]; % 三级模型中心磁场 (1,80)

    [Mutual_Inductance_3_1,Mutual_Inductance_3_2,Vertical_field_3_2] =  fun_Inductance_and_Vertical_field_gama(IC.*level_3_1_Mutual_inductance_Matrix_inner,IC.*level_3_1_Mutual_inductance_Matrix_outter,IC.*level_3_2_Mutual_inductance_Matrix_inner,IC.*level_3_2_Mutual_inductance_Matrix_outter,level_3_2_Vertical_field_inner,level_3_2_Vertical_field_outter,n_1,n_2,n_3); % lv.3 和lv.1 之间的互感；lv.3和lv.2之间的互感；lv.3 在 lv.2位置产生的垂直场
    Mutual_Inductance_3_2 = Mutual_Inductance_3_2(:,T+1:T*n_1);

    M_3_0 = [Mutual_Inductance_3_1(:,1:P_0),lv0_Mutual_inductance_Matrix(end - n_2 - n_3 + 1:end,:),Mutual_Inductance_3_1(:,P_0 + A + 1 :end)];M_0_3 = M_3_0';

    [Vertical_field_2_3,Vertical_field_3_3] = fun_Vertical_field_3(level_2_3_Vertical_field_inner,level_2_3_Vertical_field_outter,level_3_3_Vertical_field_outter_outter,level_3_3_Vertical_field_inner_inner,level_3_3_Vertical_field_outter_inner,level_3_3_Vertical_field_inner_outter,n_1,n_2,n_3);

%% 数据记录
disp('展开记录文件...')
    % 记录值
    current_record = zeros(N_4 - A + A*lv0_n,loop_number);%电流记录
    time_record = zeros(1,loop_number);%时间轴
    temperature_record = zeros(N_1,loop_number);%温度记录
    voltage_record = zeros(1,loop_number);
    field_record = zeros(1,loop_number);
    Circle_Resistance_Matrix_1_record = zeros(size(Circle_Resistance_Matrix_1,1),loop_number);
    Circle_Resistance_Matrix_2_record = zeros(size(Circle_Resistance_Matrix_2,1),loop_number);
    disp('开始循环运算...')

%% 循环计算
for loop_count = 1:loop_number

    % 循环时间
        tic
        dt = 10;
        t = t + dt;
        time_record(loop_count) = t; %记录时间轴
        now_time = datestr(now,'日期yyyy-mm-dd 时间HH:MM:SS');
        disp(['磁体时间：',num2str(t),'秒','本次循环开始时间：【',now_time,' 】']);% 显示当前时间
        prepare_time = toc;tic % 记录准备时间

    % 用于环向电阻计算的向前预估电流
    if loop_count > 2
        KI_hypothesis = 2.*KI -  current_record(:,loop_count - 2);
        I_hypothesis = fun_I_format_change_01(KI_hypothesis);  
    elseif loop_count == 2
        KI_hypothesis = 2.*KI;
        I_hypothesis = fun_I_format_change_01(KI_hypothesis);  
    else
        KI_hypothesis = KI;
        I_hypothesis = I;
    end


    % 方程分为9个模块
        % [1,1]
        % 需要先计算出垂直磁场和温度
            %垂直磁场
            J = fun_To_level_2_current(I(1:N_3),Beta_up,Beta_down);% 转化为二级电流模式，方便计算垂直磁场
            VMF_I = (sqrt(sum(X_B.*I(1:N_2)).^2 + sum(Y_B.*I(1:N_2)).^2))';%每个电流元位置的垂直场 T*A 的列向量
            VMF_J = Vertical_Field_Matrix_2_1(:,:,Beta_up+1)*J;% β系列在α面上产生的垂直场 注意这个结果是T的列向量
            VMF_K = Vertical_field_3_2(:,1:T)'*I(N_3+1:N_4);% lv.3模型在lv.1上产生的垂直磁场
            VMF = VMF_I + reshape(repmat((VMF_J+VMF_K),1,A)',T*A,1);% 进行复制之后就可以相加了
            %温度
            Heating_Power = abs(fun_Current_Heat(Circle_Resistance_Matrix_1,Radial_resistance_Matrix,I(1:N_2))); %上个计数循环的发热总量
            Temperature = fun_Temperature_change(thermal_conductivity_matrix,Heating_Power); % 上个计数循环每个节点的温度
        %构成方程
            [equation_1_1,constant_1_1,Circle_Resistance_Matrix_1,Temperature]=... %输出[1,1]位置的方程系数矩阵、常数项；α面的环向电阻和温度
            level_0_equation(t,dt,...%时间项
            I(1:N_2),...%精细模型的全部电流
            VMF,Temperature,...%垂直磁场和温度
            IC.*Single_pancake_Mutual_inductance_Matrix,Radial_resistance_Matrix,...% 精细模型的自感和径向电阻
            KI(1:(2*T-2)*A + A*lv0_n),lv0_Mutual_inductance_Matrix,lv0_p,lv0_n...
            ,current_record,loop_count); % 电流；lv0互感矩阵；lv0位置，lv0每个位置数量

            
        % [2,2]
            % 系数包括：互感和自感/dt+径向电阻+环向电阻
            % 常数包括：互感和自感乘上个时间点的电流/dt+控制电流*径向电阻
            % 常数
                canstant_2_2 =  Mutual_Inductance_Matrix_2_2*I(N_2+1:N_3)./dt+...
                                repmat(fun_Current_Control(t).*Rough_model_Radial_resistance,Beta_up + Beta_down,1);
            % 互感和自感/dt
                equation_2_2 = Mutual_Inductance_Matrix_2_2./dt;
            % 加上径向电阻
                equation_2_2 = equation_2_2 + diag(repmat(Rough_model_Radial_resistance,Beta_up + Beta_down,1));
            % 环向电阻需要计算温度和垂直场
                for k = 1:(Beta_up+Beta_down)
                    Temperature_temp = fun_level_2_Tamperature_Change(I(N_2+(k-1)*T+1:N_2+k*T),level_2_thermal_conductivity_matrix,Rough_model_Radial_resistance,Circle_Resistance_Matrix_2((k-1)*T+1:k*T),fun_Current_Control(t));
                    % 垂直磁场
                    if k > Beta_up
                        k_temp = k + 1; % 在一级模型在二级模型之间，而垂直场的计算矩阵是按照实际的计算矩阵来排列的，所以需要调整所有在一级模型之下的模型的序号，而在一级模型之上的单饼，计算其垂直场时，序号不用调整
                    end
                    VMF_temp = Vertical_Field_Matrix_2_1(:,:,k_temp)*J+...
                    Vertical_field_3_2(:,(k_temp-1)*T + 1 : k_temp*T)'*I(N_3+1:N_4);% 垂直场，注意不包括单饼的垂直自场
                    % 迭代环向电阻 %环形电阻应当是一级模型的形式参与计算 迭代用的环形电阻应当是二级模型的形式
                    Circle_Resistance_Matrix_2((k-1)*T+1:k*T) = fun_Circle_Resistance_Matrix(I_hypothesis(N_2 + (k-1)*T + 1 : N_2 + k*T),VMF_temp,Temperature_temp);
                end
                % 在方程里面加上径向电阻
                equation_2_2 = equation_2_2 + diag(Circle_Resistance_Matrix_2);
            
        % [3,3]
            equation_3_3 = IC.*level_3_Mutual_inductance_Matrix./dt +...% 一阶微分
            diag(level_3_Radial_Resistance_Matrix);% 电阻径向电阻
            canstant_3_3 = level_3_Radial_Resistance_Matrix'*fun_Current_Control(t) +... % 径向电阻乘控制电流
            (IC.*level_3_Mutual_inductance_Matrix*I(N_3 + 1 :N_4))./dt; % 一阶微分

            % 在方程里面加上径向电阻
            Circle_Resistance_Matrix_3 = fun_level_3_Circle_resistance(J,I_hypothesis(N_3+1:N_4),Vertical_field_2_3,Vertical_field_3_3,n_2,n_3);
            equation_3_3 = equation_3_3 + diag(Circle_Resistance_Matrix_3);
        % [1,2]
            %方程只有基尔霍夫电压方程的部分，也就是T*A-A行，T*(UP+DN)列，但是在构成矩阵的时候还是要有2*T*A - A行N_2
            equation_1_2 = zeros(N_2 - A + A*lv0_n,T*(Beta_up + Beta_down));
            canstant_1_2 = zeros(N_2 - A + A*lv0_n,1);
            % k = 1 第一个方程描述β面所有电流在α面中 2 ~ A+1 电流元上产生的电压 
            equation_1_2(N_1+1,:) = sum(M_0_2(2:A+1,:))./dt;
            %equation_1_2(N_1+2:N_2,:) = (Mutual_Inductance_Matrix_2_1(2:T*A-A,:) - Mutual_Inductance_Matrix_2_1(2+A:T*A,:))./dt;
            
            % k < P_0 - A 正常构成方程
            equation_1_2(N_1 + 2 : N_1+ P_0 - A,:) = (M_0_2(2:P_0 - A,:) - M_0_2(2 + A:P_0,:))./dt;
            % k < P_0 外侧的电压选择最上层的电压
            equation_1_2(N_1 + P_0 - A + 1 : N_1 + P_0,:) = (M_0_2(P_0 - A + 1:P_0,:) - M_0_2(P_0 + 1:lv0_n:P_1,:))./dt;
            % k < P_1 根据内侧的位置数量构成方程
            equation_1_2(N_1 + P_0 + 1: N_1 + P_1,:) = (M_0_2(P_0 + 1:P_1,:) - kron(M_0_2(P_1 + 1:P_1 + A,:),ones(lv0_n,1)))./dt;
            % k < P_2 - A 正常构成方程P_2
            equation_1_2(N_1 + P_1 + 1:N_1 + P_2 - A,:) = (M_0_2(P_1 + 1:P_2 - A,:) - M_0_2(P_1 + 1 + A:P_2,:))./dt;
            % 常数项
            canstant_1_2(N_1+1:end)  = equation_1_2(N_1+1:end,:)*I(N_2 + 1 : N_3);
        % [1,3]
            %方程只有基尔霍夫电压方程的部分，也就是T*A-A行，n_2 + n_3列,，但是在构成矩阵的时候还是要有N_2行
            equation_1_3 = zeros(N_2 - A + A*lv0_n,n_2 + n_3);
            canstant_1_3 = zeros(N_2 - A + A*lv0_n,1);
            % k = 1
            equation_1_3(N_1+1,:) = sum(M_0_3(2:A+1,:))./dt;
            % k < P_0 - A 正常构成方程
            equation_1_3(N_1 + 2 : N_1+ P_0 - A,:) = (M_0_3(2:P_0 - A,:) - M_0_3(2 + A:P_0,:))./dt;
            % k < P_0 外侧的电压选择最上层的电压
            equation_1_3(N_1 + P_0 - A + 1 : N_1 + P_0,:) = (M_0_3(P_0 - A + 1:P_0,:) - M_0_3(P_0 + 1:lv0_n:P_1,:))./dt;
            % k < P_1 根据内侧的位置数量构成方程
            equation_1_3(N_1 + P_0 + 1: N_1 + P_1,:) = (M_0_3(P_0 + 1:P_1,:) - kron(M_0_3(P_1 + 1:P_1 + A,:),ones(lv0_n,1)))./dt;
            % k < P_2 - A 正常构成方程P_2
            equation_1_3(N_1 + P_1 + 1:N_1 + P_2 - A,:) = (M_0_3(P_1 + 1:P_2 - A,:) - M_0_3(P_1 + 1 + A:P_2,:))./dt;
            %常数项
            canstant_1_3(N_1+1:end) = equation_1_3(N_1+1:end,:)*I(N_3+1:N_4);


        % [2,1]
            % lv.1的电流元在lv.2上产生的电压
            equation_2_1 = zeros(T*(Beta_up+Beta_down),N_2 - A + A*lv0_n);
            
            equation_2_1(:,1:N_0) = M_2_0(:,1:N_0)./dt;
            canstant_2_1 = equation_2_1*KI(1:N_2 - A + A*lv0_n);
        % [2,3]
            % lv.3的电流元在lv.2上产生的电压
            equation_2_3 = Mutual_Inductance_3_2'./dt;
            canstant_2_3 = Mutual_Inductance_3_2'*I(N_3 + 1:N_4)./dt;
        % [3,1]
            % lv.1的电流元在lv.3上产生的电压
            equation_3_1 = zeros((n_2+n_3),N_2 - A + A*lv0_n);
            
            equation_3_1(:,1:N_0) = M_3_0(:,1:N_0)./dt;
            canstant_3_1 = equation_3_1*KI(1:N_2 - A + A*lv0_n);
        % [3,2]
            % lv.2的电流元在lv.3上产生的电压
            equation_3_2 = Mutual_Inductance_3_2./dt;
            canstant_3_2 = Mutual_Inductance_3_2*I(N_2 + 1:N_3)./dt;

    % 组合成为一个方程组
        Equation = [equation_1_1,equation_1_2,equation_1_3;...
        equation_2_1,equation_2_2,equation_2_3;...
        equation_3_1,equation_3_2,equation_3_3];
        Canstant = [constant_1_1 + canstant_1_2 + canstant_1_3;...
        canstant_2_1 + canstant_2_2 + canstant_2_3;...
        canstant_3_1 + canstant_3_2 + canstant_3_3];

    % 解方程
    organize_time = toc;tic
    KI = real(Equation\Canstant);
        %% 记录  
        current_record(:,loop_count) = KI;%记录电流
        I = fun_I_format_change_01(KI); 
        temperature_record(:,loop_count) = Temperature;%记录温度
        now_voltage = level_3_Radial_Resistance_Matrix*(fun_Current_Control(t) - I(N_3+1:N_4)) + fun_Voltage_Multi_pancake(I(1:N_3),Radial_resistance_Matrix,Rough_model_Radial_resistance,t,Beta_up,Beta_down);% 总电压计算
        voltage_record(:,loop_count) = now_voltage;%记录总电压
        now_field = level_3_Central_magnetic_field_Matrix*I(N_3+1:N_4) + fun_Field_Multi_pancake(J,level_2_Central_Field_Matrix);% 总磁场计算
        field_record(:,loop_count) = now_field;%记录中心磁场
        Circle_Resistance_Matrix_1_record(:,loop_count) = Circle_Resistance_Matrix_1;
        Circle_Resistance_Matrix_2_record(:,loop_count) = Circle_Resistance_Matrix_2;
        screen_current = lv0_n*(sum(KI(P_0 + 1:lv0_n:P_1)) - sum(KI(P_0 + lv0_n:lv0_n:P_1)))/A;
    %% 显示
        solve_time = toc;tic
        disp(['准备时间',num2str(prepare_time),'秒，构成方程时间',num2str(organize_time),'秒，解方程时间',num2str(solve_time),'秒，本次循环时间',num2str(prepare_time + organize_time + solve_time),'秒']);
        disp(['预估电流精准度',num2str(sum(I_hypothesis./I)./size(I,1)),';屏蔽电流：',num2str(screen_current),'A;控制电流：',num2str(KI(1)),'A']);

        disp(['V=',num2str(now_voltage),' | B=',num2str(now_field)])
   
end
%% 保存记录
save([this_path,'/record/','record_',programdate,'.mat'],'time_record','current_record','temperature_record','voltage_record','field_record','Circle_Resistance_Matrix_1_record','Circle_Resistance_Matrix_2_record');