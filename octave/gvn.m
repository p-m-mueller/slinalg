function [y, cs, sn] = gvn(x)
  y = x;
  r = sqrt(x(1)*x(1) + x(2)*x(2));
  cs = x(1) / r;
  sn = x(2) / r;
  y(1) = r;
  y(2) = 0.0;
endfunction