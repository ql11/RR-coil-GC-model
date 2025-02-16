function varargout = Attitude(varargin)
% Attitude - 线圈基础属性，按照字典赋值
%% 使用示例：[MA,MB,MC] = Attitude('Clearance within DP','Coil number','Tape width')

%%  字典内容 myDict('key') = value; % 注释
% 初始化字典并定义
myDict = containers.Map({'Coil number'}, {2});

% 仿真基础内容
myDict('u0') = 4*pi*1e-7; %真空磁导率
myDict('N of divisions') = 8; % 每匝分为8个元素
myDict('Default time step') = 0.01;% 时间步长

% 磁体结构特征
myDict('Tape width') = 6e-3; % 带材宽度6mm
myDict('Thickness per turn') = 0.21e-3; % 每匝厚度
myDict('Pole pitch') = 480e-3; % 极距
myDict('SP N') = 406; % 单饼匝数
myDict('rated current') = 210; % 额定电流
myDict('N of DPs') = 4; % 单级双饼个数
myDict('N of coils') = 2; % 线圈个数
myDict('Edge width') = 85e-3; % 边宽
myDict('Total width') = 410e-3; % 总宽度
myDict('Total height') = 380e-3; % 总高度
myDict('Inner fillet radius') = 40e-3; % 圆角半径
myDict('Outer fillet radius') = 125e-3; % 圆角外径
myDict('Clearance within DP') = 1.4e-3; % 饼内间隙
myDict('Clearance between DP') = 3.4e-3; % 饼间间隙
myDict('Length x-axis') = 160e-3; % x轴直线段长度
myDict('Length y-axis') = 130e-3; % y轴直线段长度

% 材料、电磁特性
myDict('MMF') = 682e3; % 磁动势
myDict('Basic critical current') = 55; % 默认临界电流50A
myDict('Basic radial resistance') = 0.1; % 默认径向电阻1Ω


% 默认环境特征
myDict('Default magnetic field') = 0; % 默认磁场0T
myDict('Default temperature') = 77; % 默认温度77K

% 控制电流特征
myDict('Target current') = 4; % 目标电流
myDict('Excitation rate') = 2; % 励磁速率




%% 调用字典输出
keysToRetrieve = varargin; % 转化后格式{'name','tape width'}
valueOut = zeros(size(keysToRetrieve));

for i = 1:length(keysToRetrieve)
    keyTemp = keysToRetrieve{i};
    if isKey(myDict, keyTemp)
        valueOut(i) = myDict(keyTemp);
    else
        fprintf('Key "%s" does not exist in the attitute dictionary.\n', keyTemp);
    end
end

%% 转化输出格式
varargout = num2cell(valueOut);

end
