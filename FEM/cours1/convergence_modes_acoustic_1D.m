%
%   convergence_modes_acoustic_1D.m
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
% with a 1D linear FEM model.
%

clear all
close all
clc

% Nombre d'elements lineaires
nb_elem_max=20;
% Longueur de la cavite
L=1;
k_analytical=pi/L
% Nombre de modes calcules
nb_modes=5;


for nb_elem=1:nb_elem_max

    % Taille des elements
    size_elem=L/nb_elem;
    % Matrices elementaires
    m_elem_lin=(1/size_elem)*[1 -1;-1 1];
    k_elem_lin=(size_elem/6)*[2  1; 1 2];
    k_elem_quad=(size_elem/30)*[4 2  -1;2 16 2;-1 2 4];
    m_elem_quad=(1/(3*size_elem))*[7 -8 1;-8 16 -8;1 -8 7];


    % Creation de la matrice globale
    M_global=sparse(nb_elem+1,nb_elem+1);
    K_global=sparse(nb_elem+1,nb_elem+1);

    % Assemblage
    for ie=1:nb_elem
        dof_1=ie;
        dof_2=ie+1;
        M_global(dof_1:dof_2,dof_1:dof_2)=M_global(dof_1:dof_2,dof_1:dof_2)+m_elem_lin;
        K_global(dof_1:dof_2,dof_1:dof_2)=K_global(dof_1:dof_2,dof_1:dof_2)+k_elem_lin;
    end

    % Calcul des valeurs propres
    temp=sqrt(eigs(M_global,K_global,2,'sa'));
    k_FEM_lin(nb_elem)=max(temp);

    % Creation de la matrice globale
    M_global=sparse(2*nb_elem+1,2*nb_elem+1);
    K_global=sparse(2*nb_elem+1,2*nb_elem+1);

    % Assemblage
    for ie=1:nb_elem
        dof_1=1+2*(ie-1);
        dof_2=dof_1+2;
        M_global(dof_1:dof_2,dof_1:dof_2)=M_global(dof_1:dof_2,dof_1:dof_2)+m_elem_quad;
        K_global(dof_1:dof_2,dof_1:dof_2)=K_global(dof_1:dof_2,dof_1:dof_2)+k_elem_quad;
    end

    temp=sqrt(eigs(M_global,K_global,2,'sa'));
    k_FEM_quad(nb_elem)=max(temp);

end


figure
set(gca,'Fontsize',15)
plot(k_FEM_lin,'.-','Linewidth',3,'Markersize',25)
hold on
plot(k_FEM_quad,'r.-','Linewidth',3,'Markersize',25)
legend('Linear el.','Quadratic el.')
print -djpeg 'convergence_1'


figure
set(gca,'Fontsize',15)
plot([2:nb_elem_max+1],k_FEM_lin,'.-','Linewidth',3,'Markersize',25)
hold on
plot([3:2:(2*nb_elem_max+1)],k_FEM_quad,'r.-','Linewidth',3,'Markersize',25)
legend('Linear el.','Quadratic el.')
print -djpeg 'convergence_2'


figure
set(gca,'Fontsize',15)
loglog([2:nb_elem_max+1],k_FEM_lin-k_analytical,'-','Linewidth',3,'Markersize',25)
hold on
loglog([3:2:(2*nb_elem_max+1)],k_FEM_quad-k_analytical,'r-','Linewidth',3,'Markersize',25)
legend('Linear el.','Quadratic el.')
print -djpeg 'convergence_3'


