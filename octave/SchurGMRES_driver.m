
na = 12; n = 4;

A = -diag(ones(na-1,1),-1) + 2*diag(ones(na,1)) - diag(ones(na-1,1),1);
B = zeros(na,n);
for i = 1:n
  B(i:na,i) = 1;
endfor


% Manufactured solution
x = ones(na,1);
y = ones(n,1);

f = A * x + B * y;
g = B' * x;

% Initial guess 
x0 = zeros(na,1);
y0 = zeros(n,1);

maxIter = 3;
tol = 1.0e-16;

[x, y, iter, res] = mySchurGMRES(A, B, f, g, x, y, maxIter, tol);

errVec = [f;g] - [A, B; B', zeros(n)] * [x; y];
err = norm(errVec);

printf("Iterations: %g\n", iter);
printf("Residual:   %e\n", res);
printf("Error:      %e\n", err);