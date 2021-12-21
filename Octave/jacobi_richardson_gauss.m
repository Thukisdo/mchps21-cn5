1;

function b = set_inital_conditions(n, T0, T1)
  b = zeros(n);
  b(1) = T0;
  b(n) = T1;
endfunction

function [x, norms] = jacobi(A, x, b, epsilon)
  res_norm = 1;
  norms = [];
  while (res_norm > epsilon)
    res = (b - A * x);
    x = x + 0.5 * res;
    res_norm = norm(res, 2)  / norm(x, 2);
    norms = [norms; res_norm];
  endwhile
  
endfunction

function [x, norms] = gauss_seidel(A, x, b, epsilon)
  M = inv(tril(A))
  n = length(b);
  res_norm = 1;
  norms = [];
  while (res_norm > epsilon)
    res = (b - A * x);
    x = x + M * res;
    res_norm = norm(res, 2)  / norm(x, 2);
    norms = [norms; res_norm];
  endwhile
endfunction


n = 50;
A = zeros(n);
A(1:1+n:n*n) = 2;
A(n+1:1+n:n*n) = -1;
A(2:1+n:n*n-n) = -1;
T0 = -5;
T1 = 5;

b = set_inital_conditions(n, T0, T1);
epsilon = 1e-8;

[x, norms] = jacobi(A, zeros(n, 1), b, epsilon);
[x2, norms_gs] = gauss_seidel(A, zeros(n, 1), b, epsilon);

x_max = max([length(norms), length(norms_gs), length(norms_rs)]);

hold on
set(gca, 'YScale', 'log')
xlabel("iteration")
ylabel("Relative residual error")
grid on
plot(norms, "linestyle", "-")
plot(norms_gs, "linestyle", "-")
line([0 x_max], [epsilon epsilon])
legend("Jacobi", "Gauss seidel", "epsilon")
hold off