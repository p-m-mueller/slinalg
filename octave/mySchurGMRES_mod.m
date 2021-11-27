function [x, y, iter, res] = mySchurGMRES_mod(A, b1, b2, b3, b4, f, g, x, y, maxIter, tol)
    
  n = 4;
  na = size(A)(1);
  
  m = maxIter;
  reltol = tol * norm(g);

  % Start Schur-GMRES
  Q = zeros(n,m+1);
  H = zeros(m+1,m);
  beta = zeros(m+1,1);
  cs = zeros(m,1);
  sn = zeros(m,1);
  c = zeros(n,1);
  w = zeros(na,1);
  v = zeros(na,1);

  % Determine right hand side
  [w, iter1, res1] = myGMRES(A, f, w, na, tol);
  %b = g - B' * w;
  b = zeros(n,1);
  b(1) = g(1) - dot(b1, w);
  b(2) = g(2) - dot(b2, w);
  b(3) = g(3) - dot(b3, w);
  b(4) = g(4) - dot(b4, w);

  % Compute initial residual
  %u = B * y;
  u = zeros(na,1);
  for i = 1:na
    u(i) = b1(i) * y(1) + b2(i) * y(2) + b3(i) * y(3) + b4(i) * y(4);
  endfor

  %v = A \ u;
  [v, iter1, res1] = myGMRES(A, u, v, na, tol);
  %r = b + B' * v;
  r = zeros(n,1);
  r(1) = b(1) + dot(b1, v);
  r(2) = b(2) + dot(b2, v);
  r(3) = b(3) + dot(b3, v);
  r(4) = b(4) + dot(b4, v);

  
  beta(1) = norm(r);
  Q(:,1) = r / beta(1);
  res = abs(beta(1));

  
  for j = 1:m
    
    % Arnoldi
    %u = B * Q(:,j);
    for i = 1:na
      u(i) = b1(i) * Q(1,j) + b2(i) * Q(2,j) + b3(i) * Q(3,j) + b4(i) * Q(4,j);
    endfor
    %v = A \ u;
    v = zeros(na,1);
    [v, iter, res] = myGMRES(A, u, v, na, tol);
%    r = - B' * v;
     r(1) = -dot(b1, v);
     r(2) = -dot(b2, v);
     r(3) = -dot(b3, v);
     r(4) = -dot(b4, v);
    
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
    
    if res <= reltol;
      break;
    endif
  endfor
  iter = j;
 
  z = zeros(j,1);
  for i = j:-1:1
    z(i) = (beta(i) - dot(H(i,i:j),z(i:j))) / H(i,i);
  endfor

  y = y +  Q(:,1:j)*z;
 
  % u = f - B * y;
  for i = 1:na
    u(i) = f(i) - b1(i) * y(1) - b2(i) * y(2) - b3(i) * y(3) - b4(i) * y(4);
  endfor
  u
  [x, iter2, res2] = myGMRES(A, u, x, na, tol);

endfunction