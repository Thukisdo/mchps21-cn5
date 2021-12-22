
function res = csr_mul(A, ci, ri, x)
  
  n = size(x, 1)
  res = zeros(n, 1)
  
  for i = 1 : n
    
    s = ri(i)
    e = ri(i + 1)
    for j = s : e - 1
      c = ci(j)
      res(i) = res(i) + A(j) * x(c)
    endfor
  endfor
  
endfunction