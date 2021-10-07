
n = 5;

disp("Exact solution :");
x = ones(n,1)

A = -diag(e(1:n-1),-1) + 2*diag(e) - diag(e(1:n-1),1);
b = A*x;
x0 = zeros(n,1);

maxIter = 4;
tol = 1.0e-16;

[x, iter, res] = myGMRES(A, b, x0, maxIter, tol);

printf("Iterations: %g\n", iter);
printf("Residual:   %E\n", res); 

disp("Computed solution: ");
x