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
rep=3;
num_procesos=50; % Número de procesos simulados para cada (M,T) con los mismos parámetros (c,beta)

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
coste_obj=zeros(nM*nT,3,num_dispositivos);


for m=1:num_dispositivos
    ind=2;
    ind_coste=1;
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
                coste_reemp(n)=sum(A(ind-rep:ind-rep-1+find(A(ind-rep:end,3,m)>0,1),1,m));
                tiempo_reemp(n)=A(find(A(:,3,m)>0,1),3,m);
            end
            coste_obj(ind_coste,:,m)=[mean(coste_reemp)/mean(tiempo_reemp),M(i),T(j)];
            ind_coste=ind_coste+1;
        end
    end
end

% ---------------- GRÁFICA ----------------------



