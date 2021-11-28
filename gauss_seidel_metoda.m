function [x] = gauss_seidel_metoda(gamma)

format longG
dimension = 20;

x = zeros(dimension,1);

Ls = diag(ones(1,(dimension-1))*(-1), 1);
Us = diag(ones(1, (dimension-1))*(-1), -1);
D = diag(ones(1, dimension)*gamma);

b = ones(dimension,1)*(gamma - 2);
b(1) = b(1)+1;
b(dimension) = b(dimension)+1;
A = (Ls + Us + D);

if (~issymmetric(A))
    error("Matrix is not symetric!");
elseif (~all(eig(A) > 0)) 
    error ("Matrix is not positive definite!");
end

iterations = 0;

while (norm(A*x - b) / norm(b)) >= (10 ^-6)
    x_k = (D+Ls)\((D+Ls-A)*x + b);
    x = double(x_k);
    iterations = iterations+1;
end

display(iterations)
