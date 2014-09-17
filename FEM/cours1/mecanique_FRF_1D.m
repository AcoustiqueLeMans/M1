%
%   modes_acoustic_1D.m
%
% Copyright (C) 2014 LAUM UMR CNRS 6613 (France)
% 	Olivier DAZEL <olivier.dazel@univ-lemans.fr>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%
% This program computes the first eigenvalues of the 1D acoustic cavity
% with a 1D linear FEM model and show convergence for different modes.
%

clear all
close all
clc

% Nombre d'elements lineaires
nb_elem_max=20;
% Longueur de la cavite
L=1;
E=1;
rho=1;
% Nombre de modes calcules
nb_elem=20;


frequence_min=1;
frequence_max=5;

nb_frequences=200;

frequence=linspace(frequence_min,frequence_max,nb_frequences);


size_elem=L/nb_elem;
% Matrices elementaires
k_elem_lin=E*(1/size_elem)*[1 -1;-1 1];
m_elem_lin=rho*(size_elem/6)*[2  1; 1 2];
m_elem_quad=rho*(size_elem/30)*[4 2  -1;2 16 2;-1 2 4];
k_elem_quad=E*(1/(3*size_elem))*[7 -8 1;-8 16 -8;1 -8 7];


% Creation de la matrice globale
M_global_lin=sparse(nb_elem+1,nb_elem+1);
K_global_lin=sparse(nb_elem+1,nb_elem+1);
F_global_lin=zeros(nb_elem+1,1);
F_global_lin(1)=-1;

% Assemblage
for ie=1:nb_elem
    dof_1=ie;
    dof_2=ie+1;
    M_global_lin(dof_1:dof_2,dof_1:dof_2)=M_global_lin(dof_1:dof_2,dof_1:dof_2)+m_elem_lin;
    K_global_lin(dof_1:dof_2,dof_1:dof_2)=K_global_lin(dof_1:dof_2,dof_1:dof_2)+k_elem_lin;
end


% Creation de la matrice globale
M_global_quad=sparse(2*nb_elem+1,2*nb_elem+1);
K_global_quad=sparse(2*nb_elem+1,2*nb_elem+1);
F_global_quad=zeros(2*nb_elem+1,1);
F_global_quad(1)=-1;
% Assemblage
for ie=1:nb_elem
    dof_1=1+2*(ie-1);
    dof_2=dof_1+2;
    M_global_quad(dof_1:dof_2,dof_1:dof_2)=M_global_quad(dof_1:dof_2,dof_1:dof_2)+m_elem_quad;
    K_global_quad(dof_1:dof_2,dof_1:dof_2)=K_global_quad(dof_1:dof_2,dof_1:dof_2)+k_elem_quad;
end

for i_f=1:length(frequence)
    omega=2*pi*frequence(i_f);
    u_lin=(K_global_lin-omega^2*M_global_lin)\F_global_lin;
    norm_L2_lin(i_f)=u_lin'*M_global_lin*u_lin;
    %norm_L2_lin(i_f)=u_lin(1);
    u_quad=(K_global_quad-omega^2*M_global_quad)\F_global_quad;
    norm_L2_quad(i_f)=u_quad'*M_global_quad*u_quad;
    %norm_L2_quad(i_f)=u_quad(1);
end


figure
semilogy(frequence,abs(norm_L2_lin))
hold on
semilogy(frequence,abs(norm_L2_quad),'r')




