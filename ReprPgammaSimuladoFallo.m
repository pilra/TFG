% El fallo se detecta en el momento en el que ocurre
L=6; M=3; T=6;

rep=3;

% Costes
cc=50; % coste por mantenimiento correctivo
cp=25; % coste por mantenimiento preventivo
ci=5; % coste por inspeccion

% Parámetros proceso gamma
c=2; beta=1/4;
Tf=50; N=30;


A=zeros(1,4);
P=zeros(1,2);

for k=1:rep
    [coste,deg,tiempo,insp,p]=cost3(P(size(P,1),1),T,M,L,cc,cp,ci,c,beta,Tf,N);
    P=cat(1,P,[p(:,1) p(:,2)]);
end
P=cat(1,P,[P(size(P,1),1) 0]);


% Buscamos los puntos donde se realizan mantenimientos
mant=find(P(:,2)==0);
mant=mant(3:size(mant,1))-1;
m_corr=0; m_prev=0;
for i=1:size(mant,1)
    if (P(mant(i),2)>=L)
        m_corr(size(m_corr,2)+1)=mant(i);
    else
        m_prev(size(m_prev,2)+1)=mant(i);
    end
end
m_corr=m_corr(2:size(m_corr,2));
m_prev=m_prev(2:size(m_prev,2));



% ----------- GRÁFICA ----------------
figure(); 
h=stairs(P(:,1),P(:,2),'Color', 'black','LineWidth',1.5,'DisplayName','Deterioro'); % Proceso gamma
hold on
h(2)=line([0,P(size(P,1),1)+3],[L,L],'Color','red','LineWidth',2,'DisplayName','L');
h(3)=line([0,P(size(P,1),1)+3],[M,M],'Color','blue','LineWidth',2,'DisplayName','M');

set(gca,'XTick',0:T:ceil(P(size(P,1),1)),'YTick',[M L]);
ylim([0 L+3]); xlim([0,P(size(P,1),1)+3]);

if (size(m_prev,2) >= 1) % Mantenimiento preventivo
    h(end+1)=plot(P(m_prev,1),P(m_prev,2),'mp','LineWidth',2.5,'DisplayName','Mantenimiento preventivo'); 
end

if (size(m_corr,2) >= 1) % Mantenimiento correctivo
    h(end+1)=plot(P(m_corr,1),P(m_corr,2),'rp','LineWidth',3,'DisplayName','Mantenimiento correctivo'); 
end

legend('-DynamicLegend','Location','northwest'); legend('boxoff')
% legend(h,{'Deterioro','L','M'}, 'Location','northwest'); legend('boxoff')
xlabel('Tiempo'); ylabel('Deterioro'); title('Proceso Gamma')