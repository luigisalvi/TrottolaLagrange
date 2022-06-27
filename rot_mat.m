% @autore: Luigi Salvi 
% @data: 27/06/2022
% @contatti: luigisalvi97@libero.it

function [Rot] = rot_mat(phi,teta,psi)
Rphi = [cos(phi),sin(phi),0; -sin(phi),cos(phi),0; 0,0,1];%PRECESSIONE: rotazione attorno a Zfisso, che porta N su Xfisso
Rteta = [1,0,0; 0, cos(teta), sin(teta); 0, -sin(teta),cos(teta)];  %NUTAZIONE:rotazione attorno ad N (asse dei Nodi e coincidente con Xfisso), che porta Zfisso su Zsolidale
Rpsi = [cos(psi),sin(psi),0; -sin(psi),cos(psi),0; 0,0,1];%ROT.PROPRIA: rotazione attorno ad Z==Zs, che porta Xs a coincidere con X
[Rot] = (Rphi*Rteta*Rpsi);
end
