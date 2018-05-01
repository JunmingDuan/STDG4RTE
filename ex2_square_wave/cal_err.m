function plot_ex1(K, n);

hold on;

for K = 1:4;
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  x1 = numer1(:,1); w1 = numer1(:,2); y1 = numer1(:,3); y2 = numer1(:,4);
  plot(x1, y1, '-o');
end
plot(x1, y2, '-r');
legend('K=1', 'K=2', 'K=3', 'K=4', 'exact');

