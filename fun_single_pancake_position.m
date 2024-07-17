function h = fun_single_pancake_position(sp, dp, c)
% 输入线圈单饼序号、双饼序号和线圈序号，输出单饼垂直位置，线圈排列方式如下图所示
%   c1          c2
%       sp1
% dp 1           4
%       sp2
% dp 2           3
% dp 3           2
%                   sp2
% dp 4           1
%                   sp1


%% 线圈属性
[s,d,w] = Attitude('Clearance within DP','Clearance between DP','Tape width');


if c == 1
    hd_0 = s + d + w.*2; %一个双饼占用的高度
    hd = hd_0.*(2.5 - dp); %需计算的双饼高度
    h = hd - (sp -2).*(s + w)./2; %sp=1时增加1/2双饼厚度，sp=2时减少1/2双饼厚度
elseif c == 2
    hd_0 = s + d + w.*2; %一个双饼的高度
    hd = hd_0.*(2.5 - dp); %需计算的双饼高度
    h = -(hd - (sp -2).*(s + w)./2); % 与c=1时相反
else
    h = 0;
end
    
end


