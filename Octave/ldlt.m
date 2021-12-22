


function A = ldlt(A)
  n = size(A, 1);
  v = zeros(n, 1);
  for i = 1 : n
     for j = 1 : i-1
        v(j) = A(i, j) * A(j, j);
     endfor
     A(i, i) = A(i, i) - A(i, 1 : i - 1) * v(1 : i - 1);
     A(i + 1 : n, i) = (A(i + 1 : n, i) - A(i + 1 : n, 1 : i - 1) * v(1 : i - 1)) / A(i, i);
   endfor
endfunction