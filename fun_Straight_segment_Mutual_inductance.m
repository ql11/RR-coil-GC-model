function M = fun_Straight_segment_Mutual_inductance(l,d)
%fun_Straight_segment_Mutual_inductance - 求解平行等长直线段的互感
%
% 互感 = fun(线段长度，线段间距)

u0 = 4*pi*1e-7;

if nargin == 1 %求解自感
    M = u0.*l./(2.*pi).*(log(2.*l./5e-6)-0.75); %直线段自感
else
    M = u0.*l./(2.*pi).*(log(2.*l./d)-1); % 两直线段互感
end

