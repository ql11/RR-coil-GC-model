function M = fun_Straight_segment_Mutual_inductance(l,d)
%fun_Straight_segment_Mutual_inductance - 求解平行等长直线段的互感
%
% 互感 = fun(线段长度，线段间距)

[u0,w] = Attitude('u0','Tape width'); % u0 = 4*pi*1e-7; w带材宽度
Ld = 2.*sqrt(w/12); % 莱尔定律考虑带材宽度

if nargin == 1 %求解自感
    M = u0.*l./(2.*pi).*(log(2.*l./w)-0.75); %直线段自感
else
    M1 = u0.*l./(2.*pi).*(log(2.*l./sqrt(d.^2 + Ld.^2))-1); % 两直线段互感公式

    if real(M1) > 0 % 大于零就不管了
        M = M1;
    else % 小于0再积分计算
        fun_M = @(x,y) u0./(4*pi)./(sqrt(d.^2+(x-y).^2));
        M = integral2(fun_M,0,l,0,l);
    end
end

