
N=1000,
for ib=1:2,
if ib==1,
    beta=0.01;
end
if ib==2,
    beta=10;
end

[a,kn,Z] = NGF_d1(N,1,beta,0);

M=nnz(a)/2;
L=diag(sum(a))-a;
k=sum(a);
N=numel(k);

lambda=real(eigs(L,N-1));

M0b=0.6;
M0a=0.1;
M0c=0.6;
m=[0.02,0.01]
il=[1.05,1.05];
gc(ib)=(sum(sqrt(1./lambda))/(N+M))^(-1);
for in=2:130,

    g=gc(ib)*(m(ib)+il(ib)^(in-2));
    
    mi=1;
    norm=N+M;
par{1}=g;
par{2}=lambda;
par{3}=mi;
par{4}=norm;
setGlobalx(par);
fun = @mass0;
x0=M0a;
x=fzero(fun,x0)
Tg(in,ib)=g;
TM(in,ib)=abs(x);
Tg(1,ib)=gc(ib);
    

end
end
TM(1,:)=0;

figure
subplot(2,1,2)
TM(1,[1,2])=[0,0];
for ic=1:2,
plot(Tg(:,ic),abs(TM(:,ic))./Tg(:,ic),'LineWidth',2)
hold on
end
xlabel('$\bf{g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
ylabel('${\bf M}_1/{\bf g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
xlim([0.5,4]);
set(gca,'FontWeight','bold','FontSize',20);
subplot(2,1,1)
for ic=1:2,
plot(Tg2(:,ic),abs(TM2(:,ic))./Tg2(:,ic),'LineWidth',2)
%loglog(Tg,TM./Tg,Tg,abs(TM2)./Tg,'LineWidth',2)
hold on
end
xlabel('$\bf{g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
ylabel('${\bf M}_0/{\bf g}$','Interpreter','Latex','FontSize',26,'FontWeight','bold')
xlim([1.0,4]);
set(gca,'FontWeight','bold','FontSize',20);
