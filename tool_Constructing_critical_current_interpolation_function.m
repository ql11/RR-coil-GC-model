% tool_Constructing_critical_current_interpolation_function 输入温度值，调用data文件夹中的数据构成临界电流插值函数
% 数据来源：https://htsdb.wimbush.eu

T = 35; %温度35K
data_path = 'D:\gitlib\RR-coil-GC-model\data'; % 定义数据文件路径
width = Attitude('Tape width'); % 调用带材宽度，SI单位米

%% 导入数据
datafile = ['Shanghai Superconductor High Field Low Temperature 2G HTS ',num2str(T),' K Angle Dependence.csv']; %定义文件名
data = readtable([data_path,'\',datafile]);% 导入数据

B = data{:,2}; % 磁场值
Theta = data{:,3}; % 磁场角度
Ic = data{:,4}.*width.*100; % 临界电流，原数据单位为A/cm

%% 进行插值
% 替换NaN值
    B = replaceNaNWithAverage(B); 
    Theta = replaceNaNWithAverage(Theta); 
    Ic = replaceNaNWithAverage(Ic); 
% 插值函数生成
func = scatteredInterpolant(B,Theta,Ic,'linear','nearest');

%% 保存插值函数以供调用
eval(['Ic_',num2str(T),'=func']);
save([data_path,'\Ic_',num2str(T),'.mat'],['Ic_',num2str(T)]); % 保存为Ic_35.mat



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

