clear all
b=load('PierreAuger_Multiplex_Coauthorship/Dataset/pierreauger_multiplex.edges');
for n=1:numel(b(:,1))
A{b(n,1)}(b(n,2),b(n,3))=b(n,4);
end
A{1}(466,466)=0;
A{2}(508,508)=0;
A{1}=1.*((A{1}+A{1}'));
A{2}=1.*((A{2}+A{2}'));

ax{1}=A{2};%weighted network layer 2
ax{2}=1.*(A{2}>0); %unweighted network layer 2



M0b=0.03;
M0a=0.03;


for ic=1:2,
a=ax{ic};
k=sum(a);
kn=k+(k==0);
N=numel(k);
M=nnz(a);
L=diag(k)-a;
lambda=real(eigs(L,N-1));
[I,J,V]=find(lambda.*(lambda>10^(-5)));
lambda=real(V);
N=nnz(k);
lambda=real(eigs(L,N-1));
[I,J,V]=find(lambda.*(lambda>10^(-5)));
lambda=real(V);

M0b=2.;
M0a=0.3;
for in=1:140,
%%betti number 1 leading to M_1
    g=0.05+0.05*(in-1);
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
Tg(in,ic)=g;
TM(in,ic)=x;
%%betti number 0 leading to M_0
g=38-0.5*in;
mi=1;
par{1}=g;
par{2}=lambda;
par{3}=mi;
par{4}=norm;
setGlobalx(par);
fun = @mass;
x0=2;
x=fzero(fun,x0)
M0b=abs(x);
Tg2(in,ic)=g;
TM2(in,ic)=(x);
end
end
figure
subplot(1,2,2)
plot(Tg(:,ic),TM(:,1)./TM(:,2),'LineWidth',2)
hold on
xlabel('$\bf{g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
ylabel('${\bf M}_1^{w}/{\bf M}_1^{u}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
xlim([0,3]);
set(gca,'FontWeight','bold','FontSize',20);


subplot(1,2,1)
plot(Tg2(:,ic),(TM2(:,1))./TM2(:,2),'LineWidth',2)
hold on
xlabel('$\bf{g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
ylabel('${\bf M}_0^{w}/{\bf M}_0^{u}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
xlim([0,30]);
set(gca,'FontWeight','bold','FontSize',20);


