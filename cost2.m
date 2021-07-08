function [coste,deg,tiempo,insp,p]=cost2(tanterior,T,M,L,cc,cp,ci,a,b,Tf,N)
%Se detecta el fallo en el momento en el que ocurre
p=procgamma(a,b,Tf,N);

p(:,1)=p(:,1)+tanterior;

sL=find(p(:,2)>=L,1);
sM=find(p(:,2)>=M,1);

insp=fix(p(sM,1)/T)+1;
mprev=find(p(:,1)>=insp*T,1);

if (p(mprev,2) >= M && p(mprev,1) < p(sL,1))
    coste=ci*insp+cp;
    
    p(mprev,1)=insp*T;
    p(mprev,2)=p(mprev-1,2);
    p=p(1:mprev,:);
else 
    coste=ci*(insp-1)+cc;
    
    p=p(1:sL,:);
end


t=min(sL,mprev);
tiempo=p(t,1);
deg=p(t,2);

return
end
