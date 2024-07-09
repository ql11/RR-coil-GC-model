function M = fun_Straight_Arc_segment_Mutual_inductance(l,r,d)
% 计算直线段和1/4圆弧段之间的互感
% ![原理图.png](https://s2.loli.net/2024/07/08/iVBdvjGuKzDsmpf.png)
% l为直线段的长度；r为圆弧的半径；d为直线段与圆弧段之间的间距，当圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
% 采用矢量磁位二重积分方法计算互感，圆弧范围为0~pi/4，微分关系为dj=RIdphi 直线段积分范围为0~l
u0 = 4*pi*1e-7;

afa = @(phi,x) phi + pi/2 - atan(d/x); % 圆心与直线段上点的连线的角度 加 圆弧段角度
s = @(phi,x) r.^2 + d.^2 + x.^2 - 2.*r.*abs(d).*cos(afa(phi,x));% 余弦定理

fun_M = @(phi,x) u0./(4*pi).*r.*cos(phi)./s(phi,x);

M = integral2(fun_M,0,pi/4,0,l);

end