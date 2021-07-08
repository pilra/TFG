function X = procgammab(c,beta,tiempo,N,b)
if nargin<4
    if size(tiempo,2)>1
        N=size(tiempo,2);
    else
        N=20;
    end
    b=1;
elseif nargin <5
    b=1;        
end

if size(tiempo,2)==1
    tiempo=linspace(0,tiempo,N);
end

dt=[0, diff(tiempo.^b)];
X=0;
for i=2:size(tiempo,2)
    X(i)=X(i-1)+gamrnd(c*dt(i),beta);
end
X=horzcat(tiempo.',X.');
return
end