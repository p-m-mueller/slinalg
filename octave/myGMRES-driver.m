
n = 5;
e = ones(n,1);

A = -diag(e(1:n-1),-1) + 2*diag(e) - diag(e(1:n-1),1);
b = A*e;
x = zeros(n,1);

m = 4;
tol = 1.0e-16;
reltol = tol * norm(b);

Q = zeros(n,m+1);
H = zeros(m+1,m);
beta = zeros(m+1,1);
cs = zeros(m,1);
sn = zeros(m,1);
c = zeros(n,1);

r = b - A*x;
beta(1) = norm(r);
Q(:,1) = r / beta(1);

printf("\n%s\t\%s\n", "Iter", "Residual");
res = abs(beta(1))
printf("%d\t%f\n", 0, res);  

for j = 1:m
  
  % Arnoldi
  r = A*Q(:,j);
  for i = 1:2
    c = Q(:,1:j)'*r;
    r = r - Q(:,1:j)*c;
    H(1:j,j) = H(1:j,j) + c;
  endfor
  H(j+1,j) = norm(r);
  Q(:,j+1) = r / H(j+1,j);
  
  % Givens
  for i = 1:j-1
    tmp = cs(i) * H(i,j) + sn(i) * H(i+1,j);
    H(i+1,j) = -sn(i) * H(i,j) + cs(i) * H(i+1,j);
    H(i,j) = tmp;
  endfor  
  [H(j:j+1,j), cs(j), sn(j)] = gvn(H(j:j+1,j));

  beta(j+1) = -sn(j) * beta(j);
  beta(j) = cs(j) * beta(j);
  res = abs(beta(j+1));
  
  printf("%d\t%f\n", j, res);  
  if res <= reltol;
    break;
  endif
endfor

y = zeros(j,1);
for i = j:-1:1
  y(i) = (beta(i) - dot(H(i,i:j),y(i:j))) / H(i,i); 
endfor

x = x +  Q(:,1:j)*y;

x






