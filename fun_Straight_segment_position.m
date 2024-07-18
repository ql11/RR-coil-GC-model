function output = fun_Straight_segment_position(d,c)
%fun_Straight_segment_position - 计算直线段坐标
%
pt = Attitude('Pole pitch');%极距
[lx,ly] = Attitude('Length x-axis','Length y-axis'); %x，y轴直线段长度

if ismember(d,[1,5]) % 竖直直线
        output = 1/2.*lx.*(d - 3)./2 + (c-1.5).*pt; 
elseif ismember(d,[3,7]) % 水平直线输出y轴坐标
        output = -1/2.*ly.*(d-5)./2; % y轴坐标和线圈位置无关
else
    output = 0;
end

end