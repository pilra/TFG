L=15; nM=10; nT=10;
M=linspace(1,10,nM); % Umbral de mantenimiento preventivo
T=linspace(1,10,nT); % Tiempo entre inspecciones
rep=10;
num_procesos=2;

% Costes
cc=50; % coste por mantenimiento correctivo
cp=25; % coste por mantenimiento preventivo
ci=5; % coste por inspeccion
cd=10; % coste por inoperatividad (por unidad de tiempo)

% Parámetros proceso gamma
c=2; beta=1/4;
Tf=50; N=30;


A=zeros(1,6);
minimo=100000;
ind=1;

for j=1:nT
    for i=1:nM
        proceso=zeros((N+1)*rep,2,num_procesos);
        for n=1:num_procesos
            deg=0; t=0; aux=1;
            for k=1:rep
                [coste,deg,tiempo,insp,p]=cost(deg,t,T(j),M(i),L,cc,cp,ci,cd,c,beta,Tf,N);
                A=cat(1,A,[coste,deg,tiempo,insp,M(i),T(j)]);
                proceso(aux:(aux+size(p,1)-1),:,n)=[p(:,1),p(:,2)];
                aux=find(proceso(:,1,n)>0,1,'last');
                deg=p(size(p,1),2); t=p(size(p,1),1);
                ind=ind+1;
            end
            % Actualizamos el mínimo y nos quedamos con el
            % mejor de los procesos
            coste_real=sum(A(ind-rep:ind-1,1))-sum(A(ind-rep:ind-2,4))*ci;
            coste_medio=coste_real/t;
            if (coste_medio<minimo(1,1))
                minimo=[coste_medio,M(i),T(j),n];
                mejor_proceso=zeros((N+2)*rep,2);
                mejor_proceso(1:size(proceso,1),:)=proceso(:,:,n);
            end
        end
        
    end
end

% Filtramos la matriz para quedarnos con las filas
% correspondientes al proceso con menor coste
filas=find(A(:,5)==minimo(1,2));
P=A(filas,:);
filas=find(P(:,6)==minimo(1,3));
P=P(filas,:);
P=P((minimo(1,4)-1)*rep+1:minimo(1,4)*rep,:);

coste=P(:,1);
deg=P(:,2);
tiempo=P(:,3);
insp=P(:,4);
insp_totales=insp(end);

aux=find(mejor_proceso(:,1)>0,1,'last');
tiempo_proceso=mejor_proceso(1:aux,1);
degradacion_proceso=mejor_proceso(1:aux,2);

coste_total=sum(coste,1)-sum(insp(1:rep-1,:),1)*ci;
tiempo_total=tiempo_proceso(end,1);


% Mantenimientos realizados
mantenimiento=zeros(rep,3);
ind=1;
for i=1:rep
    if (tiempo(i)~=0)
        mantenimiento(ind,:)=[deg(i),tiempo(i),insp(i)];
        ind=ind+1;
    end
end

%coste_medio=coste_total/rep;
%tiempo_medio=tiempo_total/rep;

