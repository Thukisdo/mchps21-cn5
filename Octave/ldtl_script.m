

T = zeros(10, 1)
L = zeros(10, 1)
for i = 10 : 100
  t = 0;
  l = 0;
  for s = 1 : 5
    A = rand(i, i);
    A = tril(A) + tril(A)';
    
    tic();
    C = ldlt(tril(A));
    t = t + toc();

    tic();
    C = mylu(A);
    l = l + toc();
  endfor
  t = t / 5;
  l = l /5;

  T = [T; t];
  L = [L; l];
endfor

hold on
grid on
xlim([15 100])
# set(gca, 'YScale', 'log')
xlabel("Input size (x*x)");
ylabel("time (s)");
plot(T)
plot(L)
legend("ldlt", "lu")
hold off



A = rand(10, 10);
A = tril(A) + tril(A)'

C = ldlt(tril(A))


L = tril(C) - diag(diag(C)) + eye(10);
D = diag(diag(C));

B = L * D * L'