function I = Operation_Current(t)
% 控制电流

I0 = 0;


[I1,V] = Attitude('Target current','Excitation rate');
t1 = (I1 - I0)/V;
It = V.*t;

% 斜坡函数
    if t <= 0
        I = I0;
    elseif t < t1
        I = It;
    else
        I = I1;
    end

end