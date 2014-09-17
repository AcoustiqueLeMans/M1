%
%   convergence_modes_acoustic_1D_several_modes.m
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
% with a 1D linear FEM model and show method convergence.
%


clear all
close all
clc

% Nombre d'elements
nb_elem_min=4;
nb_elem_max=50;
% Longueur de la cavite
L=1;

% Nombre de modes calcules
nb_modes=4;
k_analytical=[1:nb_modes]*pi/L;


for nb_elem=nb_elem_min:nb_elem_max

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
    temp=sqrt(eigs(M_global,K_global,nb_modes,'sa'));
    k_FEM_lin(nb_elem,1:nb_modes)=sort(temp)';

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

    temp=sqrt(eigs(M_global,K_global,nb_modes,'sa'));
    k_FEM_quad(nb_elem,1:nb_modes)=sort(temp)';

end


figure
set(gca,'Fontsize',15)
loglog(1,1)
hold on
i_mode=2;
loglog([nb_elem_min:nb_elem_max],k_FEM_lin(nb_elem_min:nb_elem_max,i_mode)-k_analytical(i_mode-1),'-','Linewidth',3,'Markersize',25)
loglog([nb_elem_min:nb_elem_max],k_FEM_quad(nb_elem_min:nb_elem_max,i_mode)-k_analytical(i_mode-1),'r-','Linewidth',3,'Markersize',25)
i_mode=3;
loglog([nb_elem_min:nb_elem_max],k_FEM_lin(nb_elem_min:nb_elem_max,i_mode)-k_analytical(i_mode-1),'--','Linewidth',3,'Markersize',25)
loglog([nb_elem_min:nb_elem_max],k_FEM_quad(nb_elem_min:nb_elem_max,i_mode)-k_analytical(i_mode-1),'r--','Linewidth',3,'Markersize',25)
i_mode=4;
loglog([nb_elem_min:nb_elem_max],k_FEM_lin(nb_elem_min:nb_elem_max,i_mode)-k_analytical(i_mode-1),'-.','Linewidth',3,'Markersize',25)
loglog([nb_elem_min:nb_elem_max],k_FEM_quad(nb_elem_min:nb_elem_max,i_mode)-k_analytical(i_mode-1),'r-.','Linewidth',3,'Markersize',25)



xlim([nb_elem_min nb_elem_max])

print -djpeg 'convergence_several_modes'


