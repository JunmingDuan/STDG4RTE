function plot_ex1(K, n);

figure(1)
numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
x1 = numer1(:,1); w1 = numer1(:,2); y1 = numer1(:,3); y2 = numer1(:,4);
plot(x1, y1, '-o', x1, y2, '-r');

figure(2)
hold on;
TC1 = dlmread(['ex1_Nx',num2str(n),'_K',num2str(K),'_TC.dat']);
for i = 1:size(TC1, 1);
  for j = 1:length(TC1(i,:));
    if(TC1(i,j) ~= 0)
      plot(TC1(i,j)/n, (i-1)/n,'ko');
    end
  end
end
axis([0,1,0,0.1]);

