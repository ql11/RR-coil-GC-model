% 假设有一个矩阵A
A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
 
% 获取上三角矩阵
triu_A = triu(A);
 
% 获取下三角矩阵并转置
tril_A_transposed = tril(A);
tril_A_transposed = tril_A_transposed';
 
% 合并上三角和转置后的下三角矩阵形成对称矩阵
symmetric_A = triu_A + tril_A_transposed;
 
% 输出对称矩阵
disp(symmetric_A);