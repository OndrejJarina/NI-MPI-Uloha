function [x] = jacobi_metoda(gamma)

dimension = 20;

x = zeros(dimension,1);

Ls = diag(ones(1,(dimension-1))*(-1), 1);
Us = diag(ones(1, (dimension-1))*(-1), -1);
D = diag(ones(1, dimension)*gamma);

b = ones(dimension,1)*(gamma - 2);
b(1) = b(1)+1;
b(dimension) = b(dimension)+1;
A = (Ls + Us + D);
LU = Ls+Us;

rowsSum = sum(abs(LU), 1);
columnsSum = sum(abs(LU), 2);

if (all(rowsSum >= sum(abs(D))) || all(columnsSum.' >= sum(abs(D)))) 
    error("Matrix is not diagonally dominant!"); 
end

iterations = 0;

while (norm(A*x - b) / norm(b)) >= (10 ^-6)
    x_k = D\(b-(Ls+Us)*x);
    x = double(x_k);
    iterations = iterations+1;
end

% display(iterations)
% for c = 1:10
%     x_k = D\(b-(Ls+Us)*x);
%     disp(x_k)
%     x = x_k;
% end
