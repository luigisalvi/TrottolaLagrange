% @autore: Luigi Salvi 
% @data: 27/06/2022
% @contatti: luigisalvi97@libero.it
% @reference: http://www.wirgilio.it/blog/2013/03/

function movimento(tf, time, phi,teta,psi) 
%La presente funzione prende in input il tempo finae fornito anche in input
%all'ode-solver (tf), il vettore dei tempi (time) e i vettori dei valori
%assunti dai tre angoli di eulero (phi,teta,psi). In output fornisce una
%animazione in 3D della trottola in moto di rotatorio, dalla quale è
%possibile anche apprezzare il moto di precessione e nutazione dell'asse
%della trottola, che descrive nello spazio il proprio moto. 

%inizializzazione vettore dei punti della curva 
%descritta dall'asse di rotazione
curva = zeros(length(time), 3); 

%tempo linearizzato (punti di simulazione equidistanti nel tempo)
numero_punti_simul = 0.7*length(time);
tempo_simul = linspace(0, tf, numero_punti_simul);    
fig = figure('Name','Animazione');

for i = 1:numero_punti_simul
    figure(fig);
    %evita che i vari i-esimi frame si sovrappongono, 
    %cancellandoli di volta in volta
    clf 
    hold on
    view(3) %angolo di visuale 3D
    %title(sprintf('Tempo %.1g di %g', tempo_simul(i), tf))
    
    %rotazione tramite matrice di rotazione
    rotazione = rot_mat(phi(i), teta(i), psi(i));
    
    %asse di rotazione che va da 0 a 4, ma può esser allungato anche per
    %quote negative, variando asse_down
    asse_down = rotazione * [0; 0; 0]; asse_up = rotazione * [0; 0; 4];
    plot3([asse_up(1) asse_down(1)], [asse_up(2) asse_down(2)], ...
        [asse_up(3) asse_down(3)],'k', 'linewidth', 2)
    
    %superficie di una falda di cono
    [Xcy, Ycy, Zcy] = cylinder(0.1:3, 90);
    Cyl1 = rotazione * [Xcy(1,:); Ycy(1,:); Zcy(1,:)];
    Cyl2 = rotazione * [Xcy(2,:); Ycy(2,:); Zcy(2,:)];
    Cyl3 = rotazione * [Xcy(3,:); Ycy(3,:); Zcy(3,:)];

    surf([Cyl1(1,:); Cyl2(1,:); Cyl3(1,:)], ...
        [Cyl1(2,:); Cyl2(2,:); Cyl3(2,:)], ...
        [Cyl1(3,:); Cyl2(3,:); Cyl3(3,:)],"FaceAlpha",1);
    
    %luogo dei punti percorso dall'asse di rotazione
    curva(i, :) = asse_up(:);
    plot3(curva(1:i,1), curva(1:i,2), curva(1:i,3), 'r')
    axis([-3 5 -3 5 -3 5])
    xlabel('Asse x'), ylabel('Asse y'), zlabel('Asse z')
    drawnow;
end

plot3(curva(:,1), curva(:,2), curva(:,3), 'b')

end
