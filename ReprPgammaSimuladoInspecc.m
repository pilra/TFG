% El fallo se detecta en una inspección
L=6; M=3; T=6;

rep=3;

% Costes
cc=50; % coste por mantenimiento correctivo
cp=25; % coste por mantenimiento preventivo
ci=5; % coste por inspeccion
cd=10; % coste por inoperatividad (por unidad de tiempo)

% Parámetros proceso gamma
c=2; beta=1/4;
Tf=50; N=30;


A=zeros(1,4);
P=zeros(1,2);

for k=1:rep
    [coste,deg,tiempo,insp,p]=cost(P(size(P,1),2),P(size(P,1),1),T,M,L,cc,cp,ci,cd,c,beta,Tf,N);
    P=cat(1,P,[p(2:end,1) p(2:end,2)]);
end

% Buscamos los puntos donde se realizan mantenimientos
mant=find(P(:,2)==0);
mant=mant(2:size(mant,1))-1;
m_corr=0; m_prev=0; marcador=zeros(1,2);
for i=1:size(mant,1)
    if (P(mant(i),2)>=L)
        m_corr(size(m_corr,2)+1)=mant(i);
        
        mrecta=(P(mant(i),2)-L)/(P(mant(i),1)-P(find(P(:,2)>=L,1),1));
        nrecta=P(mant(i),2)-mrecta*(P(mant(i),1));
        tfallo=(L-nrecta)/mrecta;
        marcador=cat(1,marcador,[tfallo L]);
    else
        m_prev(size(m_prev,2)+1)=mant(i);
    end
end
m_corr=m_corr(2:size(m_corr,2));
m_prev=m_prev(2:size(m_prev,2));
marcador=marcador(2:size(marcador,1),:);


% ----------- GRÁFICA ----------------
figure(); 
h=stairs(P(:,1),P(:,2),'Color', 'black','LineWidth',1.5,'DisplayName','Deterioro'); %Proceso gamma
hold on
h(2)=line([0,P(size(P,1),1)+3],[L,L],'Color','red','LineWidth',2,'DisplayName','L');
h(3)=line([0,P(size(P,1),1)+3],[M,M],'Color','blue','LineWidth',2,'DisplayName','M');

set(gca,'XTick',0:T:ceil(P(size(P,1),1)),'YTick',[M L]);
ylim([0 L+3]); xlim([0,P(size(P,1),1)+3]);

if (size(m_prev,2) >= 1) % Mantenimiento preventivo
    plot(P(m_prev,1),P(m_prev,2),'mp','LineWidth',2.5,'DisplayName','Mantenimiento preventivo');
end

if (size(m_corr,2) >= 1) % Mantenimiento correctivo
    plot(P(m_corr,1),P(m_corr,2),'rp','LineWidth',3,'DisplayName','Mantenimiento preventivo');
    scatter(marcador(:,1),marcador(:,2),80,'bx','LineWidth',2,'DisplayName','Fallo del componente');
end

legend('-DynamicLegend','Location','northwest'); legend('boxoff')
xlabel('Tiempo'); ylabel('Deterioro'); title('Proceso Gamma')