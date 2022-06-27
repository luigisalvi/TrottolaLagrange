% @autore: Luigi Salvi 
% @data: 27/06/2022
% @contatti: luigisalvi97@libero.it

clc
clear
%% PARAMETRI COSTRUTTIVI/GEOMETRICI DELLA TROTTOLA
rho = 0.60; %densità legno g/cm^3
N = 50; %numero di punti di campionamento per volume
l = 1; %lato dell'cubo di campionamento
step = l/N; %distanza tra il centroide di un campione e quello contiguo
dV = step^3; %volume del singolo campione(dl*dl*dl)
s=1;k=1; 

% variabli che registrano il valore del numeratore della relazione del 
% cdm = sum(rho*coordinata-iesima)
xA=0; yA=0; zA=0; 
% paramentri del cono: se a e b sono uguali, cono a base circolare 
% => 2 dei 3 momenti principali saranno uguali tra loro
a=8; b=8; c=40;
%somma dei dA elementi compresi nella figura moltiplicati alla densità rho
massatotale = 0; 

%% COSTRUZIONE GEOMETRICA DELLA TROTTOLA   
for i = -N:N
    for j =  -N:N
        for t = -0:N %permette di ottenere una sola falda
            xi = i*step;
            yj = j*step;
            zt = t*step;
            
            %equazione di un cono a due falde parallelo all'asse z; 
            % se secondo membro pari a 1, degenera in cilindro
            if (xi^2/a + yj^2/b  < zt^2/c) 
                xA = xA + xi*rho*dV; 
                yA = yA + yj*rho*dV;
                zA = zA + zt*rho*dV;
                massatotale = massatotale + dV*rho;     
                xdentro(s) = xi;                                           
                ydentro(s) = yj; 
                zdentro(s) = zt;
                s = s+1;     
                else
                xfuori(k) = xi;
                yfuori(k) = yj;
                zfuori(k) = zt;
                k = k+1;
             end
         end 
                
     end
end

%% DEFINZIONE DEL CENTRO DI MASSA DELLA TROTTOLA 
xb = xA/massatotale;
yb = yA/massatotale;
zb = zA/massatotale; 
ForzaPeso = -massatotale*9.81;
M = ['La massa della trottola è ', num2str(massatotale), ' kg'];
disp(M);
FP =['La forza peso applicata al centro di massa della trottola è ', num2str(ForzaPeso), ' N'];
disp(FP);

%% CALCOLO DEL TENSORE DI INERZIA E DELL'ELISSOIDE DI INERZIA ASSOCIATO ALLA TROTTOLA 
global I1 I2 I3 mgl

tensore = fun_tensore_inerzia(xdentro, ydentro, zdentro, xb, yb, zb, rho, dV);

mgl = massatotale*9.81*zb;
I1 = tensore(1,1);
I2 = tensore(2,2);
I3 = tensore(3,3);

%relazioni tra alpha, beta e gamma e i momenti polari (principali) 
xr = (1/(sqrt(I1)));
yr = (1/(sqrt(I2)));
zr = (1/(sqrt(I3)));
[Xe, Ye ,Ze] = ellipsoid(xb,yb, zb , xr, yr, zr); 


%% CONDIZIONI INIZIALI ED ODE45 
%Condizioni iniziali: vari set di valori da provare
%x0 = [0;0;pi/16;0;pi/16;1];
x0 = [4;9;0.3;0.6;60;30]; %motion 1
%x0 = [4;0.2;0.3;0.7;60;30]; %motion 2
t0=0; tf=10; %range temporale di simulazione

%ODE solver: è comunque consigliato per una maggior accuratezza dei
%risultati, utilizzare ode78 oppure ode89
[t,y] = ode45(@angoli_eulero, [t0 tf], x0); 

scelta=0;
while (scelta ~= 7)
disp('1 - PLOTTING DELLA TROTTOLA ED ELLISSOIDE DI INERZIA')
disp('2 - VISUALIZZAZIONE MATRICE DI ROTAZIONE TOTALE DA SIGMA AD S ED MOMENTI ESTERNI')
disp('3 - PLOTTING DEGLI ANGOLI DI EULERO')
disp('4 - PLOTTING DELLE VELOCITA'' ANGOLARI')
disp('5 - PLOTTING ENERGIA CINETICA E CONSERVAZIONE DELL''ENERGIA MECCANICA')
disp('6 - ANIMAZIONE TROTTOLA')
disp('7 - ESCI')
scelta = input('Digitare la scelta: ');

switch scelta 

case 1 
%% FUNZIONI DI PLOTTING DI TROTTOLA ED ELLISOIDE DI INERZIA
figure
MassaC=plot3(xb,yb,zb,'or', DisplayName='Centro di massa') ;%plotta il cdm evidenziato da una x rossa
hold on 
tr = plot3(xdentro,ydentro, zdentro,'ob', DisplayName ='Masse infinitesime');
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);
%plotta  punti interni alla curva, in verde
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); title('Trottola parametrizzata');

figure
surf(Xe,Ye,Ze, 'FaceAlpha', 0.3);
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); title('Ellissoide di inerzia di una trottola');
break

case 2
%% VISUALIZZAZIONE MATRICE DI ROTAZIONE TOTALE DA SIGMA AD S ED MOMENTI ESTERNI 
syms teta phi psi MGL
Rphi = [cos(phi),sin(phi),0; -sin(phi),cos(phi),0; 0,0,1];%rotazione attorno a Zs, che porta N su Xs
Rteta = [1,0,0; 0, cos(teta), sin(teta); 0, -sin(teta),cos(teta)]; %rotazione attorno ad N (asse dei Nodi), che porta Zs su Z
Rpsi = [cos(psi),sin(psi),0; -sin(psi),cos(psi),0; 0,0,1];%rotazione attorno ad Z==Zs, che porta Xs a coincidere con X
Rotational_Matrix = (Rphi*Rteta*Rpsi); %da Sigma(fisso) a S(solidale)
disp('La matrice trasferimento delle coordinate dal sistema Sigma-fisso al sistema S-solidale è:');
disp(Rotational_Matrix);

%Proiezione dei momenti sugli assi del sistema solidale: calcolo del valore
%di M1,M2 ed M3 tramite il corss_product - sistema simbolico
g_vett = [0; 0; -1]; %vettore gravità nel riferimento fisso;
P0_O_vett = [0;0;1]; %vettore momento nel rifermento solidale;
g_solidale = (Rotational_Matrix*g_vett); %vettore gravità riportato nel riferimento solidale
momenti_valori= MGL*cross(P0_O_vett,g_solidale); 
disp('Le proiezioni dei momenti sugli assi del riferimento solidale sono:');
disp(momenti_valori);
break

case 3
%% FUNZIONI DI PLOTTING DEGLI ANGOLI DI EULERO
figure
subplot(3,1,1)
plot(t,y(:,1), 'r',DisplayName='phi')
xlabel('Tempo(s)'); ylabel('\phi (rad)'); title('ANGOLO DI PRECESSIONE');

subplot(3,1,2)
plot(t,y(:,3), 'g',DisplayName='teta')
xlabel('Tempo(s)'); ylabel('Θ (rad)'); title('ANGOLO DI NUTAZIONE');

subplot(3,1,3)
plot(t,y(:,5), 'b',DisplayName='psi')
xlabel('Tempo(s)'); ylabel('Ψ (rad)'); title('ANGOLO DI ROTAZIONE PROPRIA');
break

case 4
%% CALCOLO E PLOTTING DELLE VELOCITA' ANGOLARI
w1 = y(:,2).*sin(y(:,3)).*sin(y(:,5))+y(:,4).*cos(y(:,5));
w2 = y(:,2).*sin(y(:,3)).*cos(y(:,5))-y(:,4).*sin(y(:,5));
w3 = y(:,2).*cos(y(:,3)) + y(:,6);
w = [w1,w2,w3];

subplot(2,2,1);
plot(t,w1, DisplayName='Velocità angolare \omega1');
xlabel('Tempo(s)'); ylabel('\omega1 (rad)'); title('VELOCITA" ANGOLARE \omega1');
subplot(2,2,2);
plot(t,w2, DisplayName='Velocità angolare \omega2');
xlabel('Tempo(s)'); ylabel('\omega2 (rad)'); title('VELOCITA" ANGOLARE \omega2');
subplot(2,1,2);
plot(t,w3, DisplayName='Velocità angolare \omega3');
xlabel('Tempo(s)'); ylabel('\omega3 (rad)'); title('VELOCITA" ANGOLARE \omega3');
break

case 5
%% PLOTTING ENERGIA CINETICA E CONSERVAZIONE DELL'ENERGIA MECCANICA

%Energia cinetica
Ttot = 1/2.*I1.*(y(:,4).^2 + ((y(:,2)).^2).*(sin(y(:,3))).^2) + 1/2.*I3.*(y(:,2).*cos(y(:,3)) + y(:,6)).^2;
plot(t,Ttot);
xlabel('Tempo(s)'); ylabel('J (joule)'); title('Energia cinetica totale - forma implicita');

%Conservazione dell'energia meccanica
U=mgl*cos(y(:,3)); %Energia potenziale gravitazionale
E = Ttot+U; %Enetgia meccanica
plot(t,Ttot, 'r');
hold on
plot(t,U,'b');
hold on
plot(t,E,'g');
xlabel('Tempo [s]'); ylabel('Joule [J]'); legend('En.Cinetica','En.Potenziale', 'En.Meccanica');
break

case 6
%% MOVIMENTO LIVE TROTTOLA 
% Calcola le coordinate che i punti della trottola con determinati 
% phi,teta e psi, precedentemente calcolati tramite le equazioni 
% del sistema con ode45. Il funzionamento della presente funzione è
% interamente basato sulla matrice di rotazione.

phi = (y(:,1));
teta = (y(:,3));
psi = (y(:,5));
time = t(:,1);

movimento(tf,time,phi,teta,psi);
break

case 7
   disp('Uscita dal programma');

    otherwise 
        disp ('Opzione non valida, digitare nuovamente.')
end
end






