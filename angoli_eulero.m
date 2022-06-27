% @autore: Luigi Salvi 
% @data: 27/06/2022
% @contatti: luigisalvi97@libero.it
% @reference: Provatidis, C.G., “Revisiting the Spinning Top,” 
% International Journal of Materials and Mechanical Engineering, 
% 1 (4). 71-88. 2012.

function [ydot] = angoli_eulero(t,y) 
% Questa funzione formalizza il sistema di equazioni differenziali (6
% equazioni in 6 incognite) a partire da due sistemi di equazioni 
% (relazioni tra le velocità e gli angoli di Eulero ed relazioni sui 
% momenti di inerzia lungo gli assi solidali al rigido e la derivata del
% momento angolare K_dot) che descrivono l'intero moto della trottola in 
% funzione degli angoli di Eulero (di nutazione, rotazione propria e 
% precessione) e delle velocità angolari. L'output della funzione è 
% espresso come segue:
% y = [phi(t), phi_dot(t), teta(t), teta_dot(t), psi(t), psi_dot(t)]

global I1 I2 I3 mgl
% Metodo per ricavare x = [x1;x2] dalle originali equazioni differenziali,
% parametrizzate e trattate come combinazione di sistemi matriciali lineari
% N.B. il sistema è costruito in modo da poter verificare il moto solo in
% condizioni di simmetria della trottola (rispetto all'asse di rotazione),
% il che implica un esemplificazione nella terza componente vettoriale 
% delle equazioni di Eulero (Mz), dove grazie ad I1=I2, risulta nullo il 
% secondo addendo (I2-I1)w1w3, da cui abbiamo I3w3_dot = 0 >> w3 costante.

w3 = y(2)*cos(y(3)) + y(6);


a13 = I1*cos(y(5))*sin(y(3));
a14 = I1*cos(y(3))*sin(y(5));
a15 = -I1*sin(y(5));
a16 = -(I2-I3)*w3*sin(y(3))*cos(y(5));
a17 = (I2-I3)*w3*sin(y(5));

a23 = -I2*sin(y(5))*sin(y(3));
a24 = I2*cos(y(5))*cos(y(3));
a25 = -I2*cos(y(5));
a26 = -(I3-I1)*w3*sin(y(3))*sin(y(5));
a27 = -(I3-I1)*w3*cos(y(5));


b1 = mgl*sin(y(3))*cos(y(5)) -a13*y(2)*y(6) -a14*y(2)*y(4) ...
    -a15*y(4)*y(6) - a16*y(2) -a17*y(4);
b2 = -mgl*sin(y(3))*sin(y(5)) -a23*y(2)*y(6) -a24*y(2)*y(4) ...
    -a25*y(4)*y(6) - a26*y(2) -a27*y(4);

A = [I1*sin(y(3))*sin(y(5)), I1*cos(y(5));...
    I2*sin(y(3))*cos(y(5)), -I2*sin(y(5))];

B = [b1 ; b2];

[x] = A^-1*B;

ydot(1,1) = y(2); 
ydot(2,1) = x(1); 
ydot(3,1) = y(4); 
ydot(4,1) = x(2); 
ydot(5,1) = w3 -y(2).*cos(y(3)); 
ydot(6,1) =-x(1)*cos(y(3)) + y(2)*y(4)*sin(y(3));
end