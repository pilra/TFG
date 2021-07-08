function X = procgamma(c,beta,Tf,N)
t=linspace(0,Tf,N);
X=0;
for i=2:N
    X(i)=X(i-1)+gamrnd(c*(t(i)-t(i-1)),beta);
end
X=horzcat(t.',X.');
return
end