
% tool_Constructing_critical_current_interpolation_function 输入温度值，调用data文件夹中的数据构成临界电流插值函数
% 临界电流 = func_Ic(温度，磁场，角度)
% 临界电流单位为6mm，温度单位为K，磁场单位为T，角度单位为°，表示磁场角度与带材法相的夹角
% 数据来源：https://htsdb.wimbush.eu



T_A = [];
B_A = [];
Theta_A = [];
Ic_A = [];



for T = 15:5:45 %温度

    data_path = 'D:\gitlib\RR-coil-GC-model\data'; % 定义数据文件路径
    width = Attitude('Tape width'); % 调用带材宽度，SI单位米

    %% 导入数据
    datafile = ['Shanghai Superconductor High Field Low Temperature 2G HTS ',num2str(T),' K Angle Dependence.csv']; %定义文件名
    data = readtable([data_path,'\',datafile]);% 导入数据

    B = data{:,2}; % 磁场值
    Theta = data{:,3}; % 磁场角度
    Ic = data{:,4}.*width.*100; % 临界电流，原数据单位为A/cm

    % 替换NaN值
    B = replaceNaNWithAverage(B); 
    Theta = replaceNaNWithAverage(Theta); 
    Ic = replaceNaNWithAverage(Ic); 

    T_A = [T_A;T.*ones(size(B))];
    B_A = [B_A;B];
    Theta_A = [Theta_A;Theta];
    Ic_A = [Ic_A;Ic];

end
%% 进行插值
% 插值函数生成
fun_Ic = scatteredInterpolant(T_A,B_A,Theta_A,Ic_A,'linear','nearest');
% 保存插值函数以供调用
save([data_path,'\fun_Ic.mat'],'fun_Ic'); % 保存为Ic.mat



%% 子函数，由于原始数据中存在一些空值，在此进行处理，用平均值代替
function vec = replaceNaNWithAverage(vec)
    % 找到 NaN 的位置
    nanIndices = isnan(vec);
    
    % 初始化替换后的向量
    replacedVec = vec;
    
    for i = 1:length(vec)
        if nanIndices(i)
            if i == 1
                replacedVec(i) = vec(i + 1);
            elseif i == length(vec)
                replacedVec(i) = vec(i - 1);
            else
                replacedVec(i) = (vec(i - 1) + vec(i + 1)) / 2;
            end
        end
    end
    
    vec = replacedVec;
end

