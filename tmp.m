if ismember(Nd_j,[1,5])  % 竖直直线
    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
    d_temp = ...
        r_j - (2*r_j + lx).*... % 根据直线段位置判断位于x、y轴
        ((abs(Nd_j - Nd_i) - 1) ~= 0) +...% 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_j
        pt.*(Nc_i - Nc_j); % 竖直直线的极距加到d上
        
    M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(ly,r_i,d_temp,0);
else  % 水平直线                   
    % 圆心远离直线段时，d>0，当圆心靠近直线段时，d<0
    d_temp = ...
        r_j - (2*r_j + ly).*... % 根据直线段位置判断位于x、y轴
        ((abs(Nd_j - Nd_i) - 1) ~= 0); % 当圆弧段不位于直线段两侧时，间距d为负值，位于两侧时间距为r_j
    M(i,j) = fun_Straight_Arc_segment_Mutual_inductance(lx,r_i,d_temp,0,(Nc_i - Nc_j));% 竖直直线的极距有额外参数
end