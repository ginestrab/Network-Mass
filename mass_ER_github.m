clear all
N=1000
 

c=4.2;

x=rand(N,N);
x=x<c/N;
x=triu(x,1);
a=x+x';
M=nnz(a)/2;
L=diag(sum(a))-a;
k=sum(a);
N=numel(k);
%L=eye(N)-diag(k.^(-1))*a;
lambda=real(eigs(L,N-1));
[I,J,V]=find(lambda.*(lambda>10^(-5)));
lambda=real(V);
betti0=N-numel(V);


for ic=1:10,
    c=1+ic;
    x=rand(N,N);
x=x<c/N;
x=triu(x,1);
a=x+x';
M=nnz(a)/2;
L=diag(sum(a))-a;
k=sum(a);
N=numel(k);

lambda=real(eigs(L,N-1));
[I,J,V]=find(lambda.*(lambda>10^(-5)));
lambda=real(V);
M0b=0.03;
M0a=0.03;
for in=1:70,
    g=0.1*1.1^(in)
    mi=(M-(N-1));
    norm=N+M;
par{1}=g;
par{2}=lambda;
par{3}=mi;
par{4}=norm;
setGlobalx(par);
fun = @mass;
x0=M0a;
x=fzero(fun,x0)
M0a=abs(x);
Tg(in,ic)=g;
TM(in,ic)=x;
g=50*0.9^(in);
 g=0.1*1.1^(in);
mi=betti0;
%mi=1;
par{1}=g;
par{2}=lambda;
par{3}=mi;
par{4}=norm;
setGlobalx(par);
fun = @mass;
%x0=abs(M0b);
x0=M0b;
x=fzero(fun,x0)
M0b=abs(x);
Tg2(in,ic)=g;
TM2(in,ic)=(x);
end
end
figure
subplot(1,2,2)
for ic=1:6,
loglog(Tg(:,ic),TM(:,ic)./Tg(:,ic),'LineWidth',2)
%loglog(Tg,TM./Tg,Tg,abs(TM2)./Tg,'LineWidth',2)
hold on
end
xlabel('$\bf{g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
ylabel('${\bf M}_1/{\bf g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
xlim([0.1,100]);
ylim([0.2,0.8]);
yticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 ])
xticks([0.1 1 10 100 ])
xtickangle(0)
set(gca,'FontWeight','bold','FontSize',20);
subplot(1,2,1)
for ic=1:6,
loglog(Tg2(:,ic),TM2(:,ic)./Tg2(:,ic),'LineWidth',2)
%loglog(Tg,TM./Tg,Tg,abs(TM2)./Tg,'LineWidth',2)
hold on
end
xlabel('$\bf{g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
ylabel('${\bf M}_0/{\bf g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
xlim([0.1,100]);
xticks([0.1 1 10 100 ])
ylim([0.001,1]);
xtickangle(0)
set(gca,'FontWeight','bold','FontSize',20);
