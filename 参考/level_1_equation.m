function [Kirchhoff_Matrix,Kirchhoff_Constant,New_Circle_Resistance_Matrix_afa,Temperature] = level_1_equation(t,dt,I,VMF,Temperature,Single_pancake_Mutual_inductance_Matrix,Radial_resistance_Matrix,I_hypothesis)
    %A_equatio - 构成精细方程，也就是总体方程的左上角
    %
    % Syntax: [Kirchhoff_Matrix,Kirchhoff_Constant,New_Circle_Resistance_Matrix_afa,Temperature] = A_equatio(t,dt,I,Single_pancake_Mutual_inductance_Matrix,Radial_resistance_Matrix,self_Vertical_Magnet_field_Matrix,beta_Vertical_Magnet_field_Matrix,thermal_conductivity_matrix)

    [T,A,~,~,~] = single_pancake_parameter(0); %单饼参数
    %% [基尔霍夫电流定律](https://gitee.com/qin_lang/img/raw/master/Picgo/20200319161300.png)
        % 赋予方程组初始值 以节点计数 一共 T*A个节点       
            Kirchhoff_Current_Matrix = zeros(T*A,(2*T-1)*A);%电流系数矩阵 列数对应电流个数，环向电流`i` T*A个，径向电流`j` T*A - A 个
            Kirchhoff_Current_Constant = zeros(T*A,1);%右侧常数矩阵

        % 方程左侧
            Kirchhoff_Current_Matrix(1:T*A,1:T*A) = Kirchhoff_Current_Matrix(1:T*A,1:T*A) + diag(ones(1,T*A));
            Kirchhoff_Current_Matrix(1:T*A-1,2:T*A) = Kirchhoff_Current_Matrix(1:T*A-1,2:T*A) - diag(ones(1,T*A-1));
            Kirchhoff_Current_Matrix(1:T*A-A,T*A+1:2*T*A-A) = Kirchhoff_Current_Matrix(1:T*A-A,T*A+1:2*T*A-A) - diag(ones(1,T*A-A));
            Kirchhoff_Current_Matrix(A+1:T*A,T*A+1:2*T*A-A) = Kirchhoff_Current_Matrix(A+1:T*A,T*A+1:2*T*A-A) + diag(ones(1,T*A-A));
        % 方程右侧
            Kirchhoff_Current_Constant(T*A) =  fun_Current_Control(t);

    
    %% [基尔霍夫电压定律](https://gitee.com/qin_lang/img/raw/master/Picgo/20200319161358.png)
        % 赋予方程组初始值 以回路计数 一共 T*A - A个节点
            Kirchhoff_Voltage_Matrix = zeros((T-1)*A,(2*T-1)*A);%电压系数矩阵 列数对应电流个数 
            Kirchhoff_Voltage_Constant = zeros((T-1)*A,1);%右侧常数矩阵
        % 径向电阻和电感部分
            % 方程左侧 电压系数矩阵
            
                Kirchhoff_Voltage_Matrix(1,T*A+1) =  -1.*Radial_resistance_Matrix(1);
                Kirchhoff_Voltage_Matrix(1,2:T*A) = Kirchhoff_Voltage_Matrix(1,2:T*A) + sum(Single_pancake_Mutual_inductance_Matrix(2:A+1,2:T*A)./dt);

                
                Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+1:2*T*A-A-1) = Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+1:2*T*A-A-1)-1.*diag(Radial_resistance_Matrix(1:T*A-A-1));
                Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+2:2*T*A-A) = Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+2:2*T*A-A) + 1.*diag(Radial_resistance_Matrix(2:T*A-A));
                
                Kirchhoff_Voltage_Matrix(2:T*A-A,1:T*A) = Kirchhoff_Voltage_Matrix(2:T*A-A,1:T*A) + Single_pancake_Mutual_inductance_Matrix(2:T*A-A,1:T*A)./dt - Single_pancake_Mutual_inductance_Matrix(2+A:T*A,1:T*A)./dt;

            % 方程右侧 常数
                Inductance_Voltage_Matrix = Single_pancake_Mutual_inductance_Matrix*I(1:T*A)./dt;
                Kirchhoff_Voltage_Constant(1) = sum(Inductance_Voltage_Matrix(2:A+1));%fun_Current_Control(t);%
                Kirchhoff_Voltage_Constant(2:T*A - A) =Inductance_Voltage_Matrix(2:T*A - A) - Inductance_Voltage_Matrix(2+A:T*A);

        % 环向电流电阻部分
            %VMF_I = (sqrt(sum(X_B.*I).^2 + sum(Y_B.*I).^2))';%每个电流元位置的垂直场 T*A 的列向量
            %VMF_J = (sqrt(sum(X_B_beta.*J).^2 + sum(Y_B_beta.*J).^2))';%每个电流元位置的垂直场 T*A 的列向量
            %VMF = VMF_I + VMF_J;
            %Heating_Power = fun_Current_Heat(Circle_Resistance_Matrix,Radial_resistance_Matrix,I); %上个计数循环的发热总量
            %Temperature = fun_Temperature_change(thermal_conductivity_matrix,Heating_Power); % 上个计数循环每个节点的温度

            New_Circle_Resistance_Matrix_afa = fun_Circle_Resistance_Matrix(I_hypothesis,VMF,Temperature);%环向电阻赋值
            
            Kirchhoff_Voltage_Matrix(1,2:A+1) = Kirchhoff_Voltage_Matrix(1,2:A+1) + (New_Circle_Resistance_Matrix_afa(2:A+1))'; % 第一个方程 第
            %k = 2 ~ T*A - A
            Kirchhoff_Voltage_Matrix(2:T*A-A,2:T*A-A) = Kirchhoff_Voltage_Matrix(2:T*A-A,2:T*A-A) + diag(New_Circle_Resistance_Matrix_afa(2:T*A-A)); %第k个方程 i(k) 的系数增加
            Kirchhoff_Voltage_Matrix(2:T*A-A,2+A:T*A) = Kirchhoff_Voltage_Matrix(2:T*A-A,2+A:T*A) -diag(New_Circle_Resistance_Matrix_afa(2+A:T*A)); %第k个方程 i(k+A) 的系数减少


    %% 组合成方程组的左侧和右侧
        Kirchhoff_Matrix = [Kirchhoff_Current_Matrix;Kirchhoff_Voltage_Matrix];
        Kirchhoff_Constant = [Kirchhoff_Current_Constant;Kirchhoff_Voltage_Constant];

    
end
