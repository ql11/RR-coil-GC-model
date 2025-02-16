function M = fun_Arc_segment_Mutual_inductance(R1,R2,afa1,afa2,bta1,bta2,dx,dy,dz)
% fun_Arc_segment_Mutual_inductance - 圆弧段互感计算
% 互感 = fun(半径1，半径2，圆弧1起始角，圆弧1终末角，圆弧2起始角，圆弧2终末角，圆弧2相对于圆弧1的相对位置x,y,z)
% 默认为1/4圆弧

[u0,w] = Attitude('u0','Tape width'); % u0 = 4*pi*1e-7; w带材宽度
Ld = 2.*sqrt(w/12); % 莱尔定律考虑带材宽度

if nargin == 6 % (R1,R2,afa1,afa2,bta1,bta2)
    dx = 0;
    dy = 0;
    dz = 0;
end
if nargin == 8 % (R1,R2,afa1,afa2,bta1,bta2,dx,dy)
    dz = 0;
end
if nargin == 2 % (R1,R2)
    dx = 0;
    dy = 0;
    dz = 0;
    afa1 = 0; bta1 = 0;
    afa2 = 1/2*pi; bta2 = 1/2*pi; 
end
if nargin == 7 % (R1,R2,afa1,bta1,dx,dy,dz)
        dz = dx;dy = bta2;dx = bta1;
        bta1 = afa2;
        afa2 = 1/2*pi + afa1; bta2 = 1/2*pi + bta1;  
end


if dz == 0 && R1 == R2 && afa1 == bta1 && afa2 == bta2

    fun_M = @(afa,bta) u0./(4*pi).*(R1.*R2).*(cos(afa - bta))./Ld;
    M = integral2(fun_M,afa1,afa2,bta1,bta2);

else
    fun_distance = @(afa,bta) sqrt(...
        (R2.*cos(bta) + dx - R1.*cos(afa)).^2 ...
        + (R2.*sin(afa) + dy - R1.*sin(afa)).^2 ...
        + dz.^2 ...
        + Ld^2);
    
    fun_M = @(afa,bta) u0./(4*pi).*(R1.*R2).*(cos(afa - bta))./fun_distance(afa,bta);
    M = integral2(fun_M,afa1,afa2,bta1,bta2);
    M = abs(M);
end

end