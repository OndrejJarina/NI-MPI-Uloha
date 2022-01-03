function [x] = jacobi_metoda(gamma)
% skript sa spusti zavolanim metody s parametrom gamma, podla ktoreho 
% je vytvorena vstupna matica A

% nastavenie formatu zobrazovania cisel
format longG
dimension = 20;

% pociatocne x je rovne nulovemu vektoru 
x = zeros(dimension,1);

% vytvorenie matice A - tri zlozky (Ls - dolna trojuholnikova matica
% Us - horna trojuholnikova matica, D - diagonalna matica)
Us = diag(ones(1,(dimension-1))*(-1), 1);
Ls = diag(ones(1, (dimension-1))*(-1), -1);
D = diag(ones(1, dimension)*gamma);
A = (Ls + Us + D);

% vektor b
b = ones(dimension,1)*(gamma - 2);
b(1) = b(1)+1;
b(dimension) = b(dimension)+1;

W = eye(dimension) - D\A;

% kontrola nutnej podmienky konvergencie
if (norm(W) >= 1)
    error("Matrix does not converge")
end

% kontrola postacujucej podmienky konvergencie
% podmienkou konvergencie je, ci matica je ostro diagonalne dominantna
rowsSum = sum(abs(Ls+Us), 1);

if ~all(rowsSum < sum(abs(D))) 
    disp("Matrix is not diagonally dominant!"); 
end

iterations = 0;
% algoritmus jacobiho metody
% podmienka while = kontrola zastavovacieho kriteria
while (norm(A*x - b) / norm(b)) >= (10 ^-6)
    % vypocet iteracie neznameho vektora x
    x_k = D\(b-(Ls+Us)*x);
    x = double(x_k);
    iterations = iterations+1;
end

display(iterations)
