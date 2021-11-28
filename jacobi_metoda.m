function [x] = jacobi_metoda(gamma)
% skript sa spusti zavolanim metody s parametrom gamma, podla ktoreho 
% je vytvorena vstupna matica A

% nastavenie formatu zobrazovania cisel
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

% vektor b
b = ones(dimension,1)*(gamma - 2);
b(1) = b(1)+1;
b(dimension) = b(dimension)+1;

% kontrola podmienky konvergencie
% podmienkou konvergencie je, ci matica je ostro diagonalne dominantna
rowsSum = sum(abs(Ls+Us), 1);
columnsSum = sum(abs(Ls+Us), 2);

if (all(rowsSum >= sum(abs(D))) || all(columnsSum.' >= sum(abs(D)))) 
    error("Matrix is not diagonally dominant!"); 
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
