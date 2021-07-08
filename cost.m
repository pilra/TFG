function [coste,deg,tiempo,insp,p]=cost(deg,t,T,M,L,cc,cp,ci,cd,c,beta,Tf,N)
% Se detecta el fallo en una inspección
p=procgammab(c,beta,Tf,N);

p(:,1)=p(:,1)+t;
p(:,2)=p(:,2)+deg;

insp_previas=fix(t/T);

sM=find(p(:,2)>=M,1);
if (isempty(sM))
    insp=fix(p(size(p,1),1)/T)-insp_previas;
    coste=ci*insp;
    tiempo=0;
    deg=p(size(p,1),2);
    return
end
sL=find(p(:,2)>=L,1);
if (isempty(sL))
    sL=size(p,1)+1;
end

insp=fix(p(sM,1)/T)+1; % inspecciones necesarias para detectar que hay que hacer mantenimiento preventivo
mprev=find(p(:,1)>=insp*T,1); % índice del instante en el que se realizaría mantenimiento preventivo
if (isempty(mprev))
    insp=insp-1-insp_previas;
    coste=ci*insp;
    tiempo=0;
    deg=p(size(p,1),2);
    return
end

if (p(mprev,2) >= M && p(mprev,2) < L)
    coste=ci*insp+cp-ci*insp_previas;
    
    p(mprev,1)=insp*T;
    if (mprev>1)
        p(mprev,2)=p(mprev-1,2);
    else
        p(mprev,2)=deg;
    end
    p=p(1:mprev,:);
    p=cat(1,p,[p(size(p,1),1),0]);
else
%     mrecta=(p(sL,2)-p(sL-1,2))/(p(sL,1)-p(sL-1,1));
%     nrecta=p(sL,2)-mrecta*(p(sL,1));
%     tfallo=(L-nrecta)/mrecta;
    
    coste=ci*insp+cc+cd*(T*insp-p(sL,1))-ci*insp_previas;
    

    if (insp*T < p(size(p,1),1) && mprev>1)
        p(mprev,1)=insp*T;
        p(mprev,2)=p(mprev-1,2);
    end
    
    p=p(1:mprev,:);
    p=cat(1,p,[p(size(p,1),1),0]);
    
end

insp=insp-insp_previas;
t=size(p,1)-1;
tiempo=p(t,1);
deg=p(t,2);

return
end
