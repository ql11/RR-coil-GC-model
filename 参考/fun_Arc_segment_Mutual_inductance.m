%fun_Arc_segment_Mutual_inductance.m
function M = fun_Arc_segment_Mutual_inductance(r1,r2,h,Ply,afa1,afa2,bta1,bta2)
% 计算互感
% 输入三个变量时计算两个不同高度不同半径同轴圆环的互感、自感；输入五个变量时计算圆弧段和圆之间的互感；输入八个变量时计算两个不同高度同轴螺旋圆弧的互感
    u0 = 4*pi*1e-7;
    if nargin == 8
        
        if h == 0 && r1 == r2 && afa1 == bta1 && afa2 == bta2
            M = (u0.*r1.*abs(afa2-afa1)./(2.*pi)).*(log(2.*r1.*abs(afa2-afa1)./5e-6) - 0.75);
        else
            R1 = @(afa) r1 + (afa - afa1)/(2*pi)*Ply; % 源电流半径
            R2 = @(bta) r2 + (bta - bta1)/(2*pi)*Ply; % 目标电流半径
            fun_M = @(afa,bta) (u0./(4*pi)).*(R1(afa).*R2(bta)).*(cos(afa - bta))./sqrt(h.^2 + R1(afa).^2 + R2(bta).^2 - 2.*R1(afa).*R2(bta).*cos(bta - afa));
            M = integral2(fun_M,afa1,afa2,bta1,bta2);
        end
    elseif nargin == 3
        if h == 0 && r1 ==r2
            M = u0*r1./4;
        else
            Lowercase_k = sqrt((4.*r1.*r2)./((r1 + r2).^2 + h.^2));
            [K,E] = ellipke(Lowercase_k.^2);% 椭圆积分
            M = u0.*sqrt(r1.*r2).*((2./Lowercase_k - Lowercase_k).*K - (2./Lowercase_k).*E); %互感
        end
    elseif nargin == 6
        Lowercase_k = sqrt((4.*r1.*r2)./((r1 + r2).^2 + h.^2));
        [K,E] = ellipke(Lowercase_k.^2);% 椭圆积分
        M = u0.*sqrt(r1.*r2).*((2./Lowercase_k - Lowercase_k).*K - (2./Lowercase_k).*E).*abs(afa2-afa1)./(2*pi); %互感
    end
end