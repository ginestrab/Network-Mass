function y=mass(M)
r = getGlobalx;
g=r{1};
lambda=r{2};
mi=r{3};
norm=r{4};
y=(M-g*mi/norm)-g*sum(M./sqrt(M^2+lambda))/norm;
