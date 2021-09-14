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
num_sim=1000; % Número de procesos simulados para cada (M,T) con los mismos parámetros (c,beta)

% Costes
cc=8000; % coste por mantenimiento correctivo
cp=550; % coste por mantenimiento preventivo
ci=25; % coste por inspeccion
cd=30; % coste por inoperatividad (por unidad de tiempo)

% Parámetros proceso gamma
c=A(:,1); 
beta=A(:,2);
Tf=1; N=20;

filas_max=nM*nT*num_sim*3;

A=zeros(filas_max,6,num_dispositivos);
coste_obj=zeros(nM*nT,3,num_dispositivos);
optimo=zeros(num_dispositivos,3);

for m=1:num_dispositivos
    ind=2;
    ind_coste=1;
    for j=1:nT
        for i=1:nM
            coste_reemp=zeros(num_sim,1);
            tiempo_reemp=zeros(num_sim,1);
            for n=1:num_sim
                deg=0; t=0; w=0;
                while A(ind-1,3,m)==0
                    [coste,deg,tiempo,insp,p]=cost(deg,t,T(j),M(i),L,cc,cp,ci,cd,c(m),beta(m),Tf,N);
                    A(ind,:,m)=[coste,deg,tiempo,insp,M(i),T(j)];
                    deg=p(size(p,1),2); t=p(size(p,1),1);
                    ind=ind+1; w=w+1;
                end
                % Coste y tiempo hasta el primer reemplazamiento
                coste_reemp(n)=sum(A(ind-w:ind-1,1,m));
                tiempo_reemp(n)=A(ind-1,3,m);
                
                A(ind,:,m)=zeros(1,6); ind=ind+1;
            end
            coste_obj(ind_coste,:,m)=[mean(coste_reemp)/mean(tiempo_reemp),M(i),T(j)];
            ind_coste=ind_coste+1;
        end
    end
    [x,c_min]=min(coste_obj(:,1,m));
    optimo(m,:)=coste_obj(c_min,:,m);
end



% ---------------- GRÁFICA ----------------------
C=zeros(nM,nT,num_dispositivos);
for m=1:num_dispositivos
    C(:,:,m)=reshape(coste_obj(:,1,m),[nM,nT]);
    
    figure();
    mesh(T,M,C(:,:,m));
    xlabel('T'); ylabel('M (mm)');
    title(['Dispositivo ' num2str(m)]);
end


