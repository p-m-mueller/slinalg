function [x, y, iter, res] = mySchurGMRES(A, B, f, g, x, y, maxIter, tol)
    
  n = size(B)(2);
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
  b = g - B' * w;

  % Compute initial residual
  u = B * y;
  %v = A \ u;
  [v, iter1, res1] = myGMRES(A, u, v, na, tol);
  r = b + B' * v;

  beta(1) = norm(r);
  Q(:,1) = r / beta(1);
  res = abs(beta(1));

  for j = 1:m
    
    % Arnoldi
    u = B * Q(:,j);
    %v = A \ u;
    v = zeros(na,1);
    [v, iter, res] = myGMRES(A, u, v, na, tol);
    r = - B' * v;
    
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

  [x, iter2, res2] = myGMRES(A, f-B*y, x, na, tol);
    
endfunction