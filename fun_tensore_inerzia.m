% @autore: Luigi Salvi 
% @data: 27/06/2022
% @contatti: luigisalvi97@libero.it

function [tensore] = fun_tensore_inerzia(xdentro, ydentro, zdentro, ...
    xb, yb, zb, rho, dV)
% La funzione calcola i momenti di inerzia principali (rispetto agli assi)
% e i momenti di inerzia secondari (centrifughi). Infine mostra a video
% la matrice simmetrica associata al tensore di inerzia, composta dei
% singoli momenti calcolati in fase di elaborazione.

Ixx = 0; % momento polare rispetto all'asse x - sigma11
for i = 1: length(xdentro) % == length (ydentro)
     %distanza del singolo dV dal cdm rispetto a z
    dr = sqrt((zdentro(i) -zb)^2 + (ydentro(i) - yb)^2);
    Ixx = dr^2 * rho * dV + Ixx;
end

Iyy = 0; % momento polare rispetto all'asse y - sigma22
for j = 1: length(xdentro) % == length (ydentro)
    %distanza del singolo dV dal cdm rispetto a y
    dr = sqrt((xdentro(j) -xb)^2 + (zdentro(j) - zb)^2); 
    Iyy = dr^2 * rho * dV + Iyy;
end

Izz = 0; % momento polare rispetto all'asse z - sigma33
for k = 1: length(xdentro) % == length (ydentro)
    %distanza del singolo dV dal cdm rispetto a z
    dr = sqrt((xdentro(k) -xb)^2 + (ydentro(k) - yb)^2); 
    Izz = dr^2 * rho * dV + Izz;
end

%momento centrifugo rispetto al piano individuato da xy -sigma12==sigma21
Ixy = 0; 
for i = 1: length(xdentro)
    rxry = ydentro(i)*xdentro(i) ;
    Ixy = - rxry * rho * dV + Ixy;
end 

%momento centrifugo rispetto al piano individuato da xz -sigma13==sigma31
Ixz = 0;
for i = 1: length(xdentro)
    rxrz = zdentro(i)*xdentro(i) ;
    Ixz = - rxrz * rho * dV + Ixz;
end

%momento centrifugo rispetto al piano individuato da yz -sigma23==sigma32
Iyz = 0;
for i = 1: length(xdentro)
    ryrz = zdentro(i)*ydentro(i) ;
    Iyz = - ryrz * rho * dV + Iyz;
end


disp('Momenti polari ottenuti sono:');
Ixx  %#ok<*NOPRT>                                                                        
Iyy
Izz

if Ixy && Ixz && Iyz < 10^-5
    ixy = round(Ixy);
    ixz = round(Ixz);
    iyz = round(Iyz);
    disp('I momenti centrifughi sono:')
    ixy 
    ixz
    iyz
else
    disp('I momenti centrifughi sono:')
    Ixy
    Ixz
    Iyz
end

tensore = [Ixx,Ixy,Ixz;Ixy,Iyy,Iyz;Ixz,Iyz,Izz];
disp('La matrice associata al tensore di inerzia è:')
disp(tensore);

end

