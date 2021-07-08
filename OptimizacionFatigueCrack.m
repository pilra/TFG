% Leemos el archivo con los mejores parámetros para cada dispositivo
fileID=fopen('fc_mejor.txt','r');
A=fscanf(fileID,'%f %f',[2 10])';
fclose(fileID);
    % Cada fila (c  beta) corresponde a un dispositivo

L=0.4; % El umbral de mantenimiento correctivo es conocido
nM=10; nT=20; 
num_dispositivos=size(A,1); % Número de dispositivos a optimizar
M=linspace(0.01,L-0.01,nM); % Umbral de mantenimiento preventivo
T=linspace(0.01,0.9,nT); % Tiempo entre inspecciones
rep=7;
num_procesos=3; % Número de procesos simulados para cada (M,T) con los mismos parámetros (c,beta)

% Costes
cc=50; % coste por mantenimiento correctivo
cp=25; % coste por mantenimiento preventivo
ci=5; % coste por inspeccion
cd=10; % coste por inoperatividad (por unidad de tiempo)

% Parámetros proceso gamma
c=A(:,1); 
beta=A(:,2);
Tf=1; N=10;

filas_max=1+nM*nT*rep;

A=zeros(filas_max,6,num_dispositivos);
minimo=repmat(100000,[num_dispositivos,4]);
mejor_proceso=zeros((N+2)*rep,2,num_dispositivos);


for m=1:num_dispositivos
    ind=2;
    
    for j=1:nT
        for i=1:nM
            proceso=zeros((N+1)*rep,2,num_procesos);
            coste_reemp=zeros(num_procesos,1);
            tiempo_reemp=zeros(num_procesos,1);
            for n=1:num_procesos
                deg=0; t=0; aux=1;
                for k=1:rep
                    [coste,deg,tiempo,insp,p]=cost(deg,t,T(j),M(i),L,cc,cp,ci,cd,c(m),beta(m),Tf,N);
                    A(ind,:,m)=[coste,deg,tiempo,insp,M(i),T(j)];
                    proceso(aux:(aux+size(p,1)-1),:,n)=[p(:,1),p(:,2)];
                    aux=find(proceso(:,1,n)>0,1,'last');
                    deg=p(size(p,1),2); t=p(size(p,1),1);
                    ind=ind+1;
                end
                % Coste hasta el primer reemplazamiento
                coste_reemp(n)=sum(A(1:find(A(:,3,m)>0,1),1,m));
                tiempo_reemp(n)=A(find(A(:,3,m)>0,1),3,m);
                
%                 % Actualizamos el mínimo y nos quedamos con el
%                 % mejor de los procesos
%                 coste_real=sum(A(ind-rep:ind-1,1,m))-sum(A(ind-rep:ind-2,4,m))*ci;
%                 coste_medio=coste_real/t;
%                 if (coste_medio<minimo(m,1))
%                     minimo(m,:)=[coste_medio,M(m,i),T(j),n];
%                     mejor_proceso(:,:,m)=zeros((N+2)*rep,2);
%                     mejor_proceso(1:size(proceso,1),:,m)=proceso(:,:,n);
%                 end
            end
            coste_obj=[mean(coste_reemp)/mean(tiempo_reemp),M(i),T(j)];
        end
    end
end


% Filtramos la matriz para quedarnos con las filas
% correspondientes al proceso con menor coste
P_aux=zeros(rep*num_procesos,6,num_dispositivos);
for m=1:num_dispositivos
    filas=find(A(:,5,m)==minimo(m,2));
    aux=A(filas,:,m);
    filas=find(aux(:,6)==minimo(m,3));
    P_aux(:,:,m)=aux(filas,:);
end

P=zeros(rep,6,num_dispositivos);
for m=1:num_dispositivos
    P(:,:,m)=P_aux((minimo(m,4)-1)*rep+1:minimo(m,4)*rep,:,m);
end
    
coste=reshape(P(:,1,:),rep,num_dispositivos);
deg=reshape(P(:,2,:),rep,num_dispositivos);
tiempo=reshape(P(:,3,:),rep,num_dispositivos);
insp=reshape(P(:,4,:),rep,num_dispositivos);
insp_totales=insp(size(insp,1),:);

tiempo_proceso=reshape(mejor_proceso(:,1,:),(N+2)*rep,num_dispositivos);
degradacion_proceso=reshape(mejor_proceso(:,2,:),(N+2)*rep,num_dispositivos);

coste_total=sum(coste,1)-sum(insp(1:rep-1,:),1)*ci;
for m=1:num_dispositivos
    aux=find(tiempo_proceso(:,m)>0,1,'last');
    tiempo_total(m)=tiempo_proceso(aux,m);
end

% Mantenimientos realizados
mantenimiento=zeros(rep,3,num_dispositivos);
for m=1:num_dispositivos
    ind=1;
    for i=1:rep
        if (tiempo(i,m)~=0)
            mantenimiento(ind,:,m)=[deg(i,m),tiempo(i,m),insp(i,m)];
            ind=ind+1;
        end
    end
end



%tiempo_medio=tiempo_total/rep;




% ---------------- GRÁFICA ----------------------
for m=1:num_dispositivos
    % Eliminamos las filas nulas del proceso
    p=find(mejor_proceso(:,1,m)>0,1,'last');
    p=mejor_proceso(1:p,:,m);
    
    figure(); 
    h=stairs(p(:,1),p(:,2),'Color', 'black','LineWidth',1.5,'DisplayName','Deterioro'); % Proceso gamma
    hold on
    h(2)=line([0,p(size(p,1),1)+0.3],[L,L],'Color','red','LineWidth',2,'DisplayName','L');
    h(3)=line([0,p(size(p,1),1)+0.3],[P(1,5,m),P(1,5,m)],'Color','blue','LineWidth',2,'DisplayName','M');

    set(gca,'XTick',0:P(1,6,m):ceil(p(size(p,1),1)),'YTick',[P(1,5,m) L]);
    ylim([0 L+0.5]); xlim([0,p(size(p,1),1)+0.3]);
    
    m_prev=0; m_corr=0;
    for i=1:size(P(:,:,m),1)
        if (P(i,3,m)>0)
            if (P(i,2,m)>P(1,5,m) && P(i,2,m)<L) % Mantenimiento preventivo
                m_prev(size(m_prev,2)+1)=i;
            elseif (P(i,2,m)>=L) % Mantenimiento correctivo
                m_corr(size(m_corr,2)+1)=i;
            end
        end
    end
    m_corr=m_corr(2:size(m_corr,2));
    m_prev=m_prev(2:size(m_prev,2));

    if (size(m_prev,2) >= 1) % Mantenimiento preventivo
        plot(P(m_prev,3,m),P(m_prev,2,m),'mp','LineWidth',2.5,'DisplayName','Mantenimiento preventivo');
    end

    if (size(m_corr,2) >= 1) % Mantenimiento correctivo
        plot(P(m_corr,3,m),P(m_corr,2,m),'rp','LineWidth',3,'DisplayName','Mantenimiento correctivo');
    end
    
    text=['T=' num2str(P(1,6,m))];
    text2=['Coste=' num2str(coste_total(m))];
    annotation('textbox',[0.7 0.7 0.1 0.18],'String',{text,text2},'FitBoxToText','on','LineStyle','none')
    legend('-DynamicLegend','Location','northwest'); legend('boxoff')
    xlabel('Cientos de miles de ciclos'); ylabel('Deterioro acumulado (mm)');
    title(['Dispositivo ' num2str(m)]);

end
