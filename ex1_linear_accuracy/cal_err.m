function plot_ex1(K, n);

format long;
hold on;

numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
x1 = numer1(:,1); w1 = numer1(:,2); y1 = numer1(:,3); y2 = numer1(:,4);
plot(x1, y1, '-o', x1, y2, '-r');

%N = [4 8 16 24 32 40 48 56 64 72 80];
%N = [10 20 40 80 160 320 640 1280]; % 2560 5120];
N = [10 20 40 80 160 320 640];
m = length(N);
err = zeros(m, 3);
min_y = zeros(m, 1);
for i = 1:m;
  n = N(i);
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  w1 = numer1(:,2); y1 = numer1(:,3); y2 = numer1(:,4);
  err(i, 1) = sqrt(sum((y1-y2).^2.*w1)/n);
  err(i, 2) = sum(abs(y1-y2).*w1)/n;
  err(i, 3) = max(abs(y1-y2));
  min_y(i) = min(y1);
end
order = zeros(m, 3);
for i = 1:m-1;
  order(i+1,:) = -log(err(i,:)./err(i+1,:))./log(N(i)/N(i+1));
end

diary table1.dat
diary on;
for j = 1:m;
  fprintf('%3d ', N(j));
  for i = 1:3
    fprintf('& %.3e & %.2f ', err(j,i), order(j,i));
  end
  fprintf('& %.3e\n', min_y(j));
  fprintf('\\\\ \n');
end
diary off;


