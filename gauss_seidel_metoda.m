function [x] = gauss_seidel_metoda(gamma)
% skript sa spusti zavolanim metody s parametrom gamma, podla ktoreho 
% je vytvorena vstupna matica A

format longG
dimension = 20;

% pociatocne x je rovne nulovemu vektoru
x = zeros(dimension,1);

% vytvorenie matice A - tri zlozky (Ls - horna trojuholnikova matica
% Us - dolna trojuholnikova matica, D - diagonalna matica)
Ls = diag(ones(1,(dimension-1))*(-1), 1);
Us = diag(ones(1, (dimension-1))*(-1), -1);
D = diag(ones(1, dimension)*gamma);
A = (Ls + Us + D);

% vytvorenie vektora b
b = ones(dimension,1)*(gamma - 2);
b(1) = b(1)+1;
b(dimension) = b(dimension)+1;

% kontrola podmienky konvergencie
% matica A musi byt symetricka a pozitivne definitna
if (~issymmetric(A))
    error("Matrix is not symetric!");
    
% vlastne cisla su vacsie ako 0
elseif (~all(eig(A) > 0)) 
    error ("Matrix is not positive definite!");
end

iterations = 0;

%algoritmus gauss-seidelovej metody
while (norm(A*x - b) / norm(b)) >= (10 ^-6)
    % vypocet iteracie neznameho vektora x
    x_k = (D+Ls)\((D+Ls-A)*x + b);
    x = double(x_k);
    iterations = iterations+1;
end

display(iterations)
