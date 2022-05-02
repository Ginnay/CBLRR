function [n_space,eigenvalues] = cal_eigenvalues(W)
%input: similarity matrix W
WB = W;
n = size(W,1);
D = diag(WB*ones(n,1));
Prw = eye(size(W)) - D^(-1/2)*WB*D^(-1/2);
if n>=1000
    No_eigs = 100;
    all_eigs = real(eigs(Prw,No_eigs,'sm'));
else
    all_eigs = real(eig(Prw));
end
    
zz = sort(abs(real(all_eigs)));  

A=1:30;
B=zz(1:30);
x=A;
y=B;
for i=1:length(x)-1
 z(i) = (y(i+1)-y(i))/(x(i+1)-x(i));
%z1(i) = y(i)/x(i);
end
n_space = find(z<0.01,2); %camp时候输出2和7

display('Number of cluster based on eigenvalues ');
display([n_space]);

eigenvalues = zz;

scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'filled');
box on;
set(gca,'LineWidth',1.5);
xlabel('i');
ylabel('Eigenvalue of graph Laplacian \lambda_i');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
%print(['Results\' folder '\EigenGap'],'-dpdf','-r300'); 
end

