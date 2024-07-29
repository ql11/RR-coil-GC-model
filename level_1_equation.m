function [Kirchhoff_Matrix,Kirchhoff_Constant] = level_1_equation(Mutual_inductance_Matrix,Radial_resistance_Array,t,dt,Magnetic_Field,Temperature,I_hypothesis)
    % level_1_equation- 构成精细方程，也就是总体方程的左上角
    %
    % Syntax: [方程系数,方程右侧常数] = 函数(互感矩阵,径向电阻数组,时间,步长,磁场数组,温度数组,预测电流)
    % level_1_equation - 构成精细方程，也就是总体方程的左上角

    %% 初始化参数
    [A,T] = Attitude('N of divisions','SP N'); %单饼参数
    if nargin == 2
        t = 10;
        
        [dt,Magnetic_Field,Temperature] = Attitude('Default time step','Default magnetic field','Default temperature');

        I_hypothesis = zeros(size(Radial_resistance_Array));
    end


    %% [基尔霍夫电流定律](https://s21.ax1x.com/2024/07/29/pkLQzNT.png)
        % 赋予方程组初始值 以节点计数 一共 T*A个节点       
            Kirchhoff_Current_Matrix = zeros(T*A,(2*T-1)*A);%电流系数矩阵 列数对应电流个数，环向电流`i` T*A个，径向电流`j` T*A - A 个
            Kirchhoff_Current_Constant = zeros(T*A,1);%右侧常数矩阵

        % 方程左侧
            Kirchhoff_Current_Matrix(1:T*A,1:T*A) = Kirchhoff_Current_Matrix(1:T*A,1:T*A) + diag(ones(1,T*A));
            Kirchhoff_Current_Matrix(1:T*A-1,2:T*A) = Kirchhoff_Current_Matrix(1:T*A-1,2:T*A) - diag(ones(1,T*A-1));
            Kirchhoff_Current_Matrix(1:T*A-A,T*A+1:2*T*A-A) = Kirchhoff_Current_Matrix(1:T*A-A,T*A+1:2*T*A-A) - diag(ones(1,T*A-A));
            Kirchhoff_Current_Matrix(A+1:T*A,T*A+1:2*T*A-A) = Kirchhoff_Current_Matrix(A+1:T*A,T*A+1:2*T*A-A) + diag(ones(1,T*A-A));
        % 方程右侧
            Kirchhoff_Current_Constant(T*A) =  Operation_Current(t);

    
    %% [基尔霍夫电压定律](https://s21.ax1x.com/2024/07/29/pkLlS4U.png)
        % 赋予方程组初始值 以回路计数 一共 T*A - A个节点
            Kirchhoff_Voltage_Matrix = zeros((T-1)*A,(2*T-1)*A);%电压系数矩阵 列数对应电流个数 
            Kirchhoff_Voltage_Constant = zeros((T-1)*A,1);%右侧常数矩阵
        % 径向电阻和电感部分
            % 方程左侧 电压系数矩阵
            
                Kirchhoff_Voltage_Matrix(1,T*A+1) =  -1.*Radial_resistance_Array(1);
                Kirchhoff_Voltage_Matrix(1,2:T*A) = Kirchhoff_Voltage_Matrix(1,2:T*A) + sum(Mutual_inductance_Matrix(2:A+1,2:T*A)./dt);

                
                Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+1:2*T*A-A-1) = Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+1:2*T*A-A-1)-1.*diag(Radial_resistance_Array(1:T*A-A-1));
                Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+2:2*T*A-A) = Kirchhoff_Voltage_Matrix(2:T*A-A,T*A+2:2*T*A-A) + 1.*diag(Radial_resistance_Array(2:T*A-A));
                
                Kirchhoff_Voltage_Matrix(2:T*A-A,1:T*A) = Kirchhoff_Voltage_Matrix(2:T*A-A,1:T*A) + Mutual_inductance_Matrix(2:T*A-A,1:T*A)./dt - Mutual_inductance_Matrix(2+A:T*A,1:T*A)./dt;

            % 方程右侧 常数
                Inductance_Voltage_Matrix = Mutual_inductance_Matrix*I(1:T*A)./dt;
                Kirchhoff_Voltage_Constant(1) = sum(Inductance_Voltage_Matrix(2:A+1));%fun_Current_Control(t);%
                Kirchhoff_Voltage_Constant(2:T*A - A) =Inductance_Voltage_Matrix(2:T*A - A) - Inductance_Voltage_Matrix(2+A:T*A);

        % 环向电流电阻部分

            New_Circle_Resistance = fun_Circle_Resistance(I_hypothesis,Magnetic_Field,Temperature);%环向电阻赋值
            
            Kirchhoff_Voltage_Matrix(1,2:A+1) = Kirchhoff_Voltage_Matrix(1,2:A+1) + (New_Circle_Resistance(2:A+1))'; % 第一个方程 第
            %k = 2 ~ T*A - A
            Kirchhoff_Voltage_Matrix(2:T*A-A,2:T*A-A) = Kirchhoff_Voltage_Matrix(2:T*A-A,2:T*A-A) + diag(New_Circle_Resistance(2:T*A-A)); %第k个方程 i(k) 的系数增加
            Kirchhoff_Voltage_Matrix(2:T*A-A,2+A:T*A) = Kirchhoff_Voltage_Matrix(2:T*A-A,2+A:T*A) -diag(New_Circle_Resistance(2+A:T*A)); %第k个方程 i(k+A) 的系数减少


    %% 组合成方程组的左侧和右侧
        Kirchhoff_Matrix = [Kirchhoff_Current_Matrix;Kirchhoff_Voltage_Matrix];
        Kirchhoff_Constant = [Kirchhoff_Current_Constant;Kirchhoff_Voltage_Constant];

    
end
