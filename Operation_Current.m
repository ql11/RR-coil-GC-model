function I = Operation_Current(t)
% 控制电流

[I1,V] = Attitude('Target current','Excitation rate');

% 斜坡函数
I = max(min(V.*t,I1),0);

end