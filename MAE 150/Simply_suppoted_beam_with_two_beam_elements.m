%This code uses two beam elements to solve a simply supported beam with
%concentrated and distributed loading. Symbolic variables are used so that
%the solution process may be more easily followed.

clear all; clc; close all;

%Define load and beam variables and parameters as symbolic
syms a E I
L = [a, 2*a]; %lengths of each beam element
a0 = 1; %in units of a (used for plotting purpurses)
L0 = eval(subs(L,a,a0));
span0 = sum(L0); %total length in units of a (used for plotting purpurses)

%Define problem size
num_elements = length(L); %number of elements
num_nodes = num_elements + 1; %number of nodes
DOF = 2*num_nodes; %degrees of freedom (two per node)

%Define loading
syms P q
necF = [2]; %nodes where external concentrated forces are applied
ecF = [P]; %corresponding magnitude of external concentrated forces
necM = []; %nodes where external concentrated moments are applied
ecM = []; %corresponding magnitude of external concentrated moments
edF = [1, 2]; %elements over which distributed forces are applied
dF = [2*q, q]; %corresponding magnitude of distributed forces

%Define reactions of supports
nYr = [1 3]; %nodes at which there is a vertical reaction
nMr = []; %nodes at which there is a moment reaction

% Visualize the undeformed structure
undeformedObj = visualize_undeformed_structure(span0);

%Visualize supports
visualize_supports(span0, a0);

%Visualize distributed loading
visualize_distributed_loading(L0, span0, q, edF, dF)

%Visualize concentrated loading
visualize_concentrated_loading(L0, necF)

%Assemble the element stiffness matrices
K = assemble_global_stiffness_matrix(num_elements, DOF, L, E, I);

%Apply external loading
F = apply_concentrated_loading(DOF, nYr, necF, nMr, necM, ecF, ecM); %Create global vector of nodal forces (and moments) accounting for concentrated forces (and moments)
B = apply_distributed_load(edF, dF, L); %Create global vector of nodal forces accounting for distributed loading

%Extract from larger system the part involving unkown generalized displacements only
[Kr, Fr, Br, iR, iud] = part_system(DOF, nYr, nMr, K, F, B);
%Kr, Fr, and Br are the reduced stiffness matrix and load vectors
%iR and iud are indices corresponding to unkown reactions and displacements in F and u respectively

%Invert reduced matrix and solve for unknown displacements
ur = inv(Kr)*(Fr+Br);

%Construct full global displacement vector
u = assemble_global_displacement_vector(DOF, iR, iud, ur);

%Solve for the unknown reactions
F = solve_for_reactions(K, u, B, nYr, nMr);

%Specifying specific values of parameters to be able to plot
P0 = -10e3; %N
q0 = -5e3/sum(subs(L,a,a0)); %N/m
E0 = 200e9; %Pa
I0 = 1/12*(0.1)^4; %m^4
u0 = eval(subs(u, [P, q, a, E, I], [P0, q0, a0, E0, I0]));

%Interpolating polynomial and visualizing deformed structure
visualize_deformation(num_elements, L0, u0, undeformedObj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%&&%%% Local functions %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function undeformedObj = visualize_undeformed_structure(span0)
figure
set(0,'defaultAxesFontName', 'Times')
undeformedObj = line([0 span0],[0 0],'color','k','linewidth',2);
axis equal
grid on
title('Undeformed and deformed beam') 
set(gca,'fontsize',14);
xlabel('{\itz} (m)');
ylabel('{\ity} (m)');
hold on
end

function visualize_supports(span0, a0)
plot(0,-0.1,'k^','MarkerSize',10)
plot(span0*a0,-0.07,'ko','MarkerSize',10)
axis([-0.5*a0, span0+0.5, -1, 1])
end

function visualize_distributed_loading(L0, span0, q, edF, dF)
N = 25; %max number of arrows if distributed load to span entire beam
for e = edF
    s = eval(dF/q); %scaling
    z1 = sum(L0(1:e-1)); %z-coord of left end of beam element
    z2 = sum(L0(1:e)); %z-coord of right end of beam element
    Ne = round((z2-z1)/span0*N); %number of arrows over element
    d = (z2-z1)/(Ne-1); %actual separation between arrows 
    for i = 1:Ne
        quiver(z1+(i-1)*d, s(e)*0.2, 0, -s(e)*0.2, 'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 5); % Disable scaling
    end
end
end

function visualize_concentrated_loading(L0, necF)
for n = necF 
    p1 = [sum(L0(1:n-1)), 0.5]; 
    p2 = [sum(L0(1:n-1)), 0];
    dp = p2 - p1;
    quiver(p1(1), p1(2), dp(1), dp(2), 'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 5); % Disable scaling
end
end

function K = assemble_global_stiffness_matrix(num_elements, DOF, L, E, I)
K = sym(zeros(DOF)); %Global stiffness matrix preallocation
for i=1:num_elements
    Keg = sym(zeros(DOF)); %Preallocation of element stiffness matrix expanded to global size
    fprintf('Individual stiffness matrix for element %d\n', i)
    Ke = E*I/L(i)^3*[ 12,        6*L(i),     -12,        6*L(i);
                      6*L(i),    4*L(i)^2,   -6*L(i),    2*L(i)^2;
                     -12,       -6*L(i),      12,       -6*L(i);
                      6*L(i),    2*L(i)^2,   -6*L(i),    4*L(i)^2];
    Keg([2*i-1:2*i+2],[2*i-1:2*i+2]) = Ke
    K = K + Keg;
end
fprintf('Global stiffness matrix for entire system\n')
K
end

function F = apply_concentrated_loading(DOF, nYr, necF, nMr, necM, ecF, ecM)
F = sym('F', [1, DOF]).'; %Define generalized forces vector as symbolic: [F1, F2=, ..., F6], i.e., [V1, M1, V2, M2, V3, M3]
izl = setdiff(1:DOF, [2*nYr-1, 2*necF-1, 2*nMr, 2*necM]); %Indices of zero load components: F2=M1=0, F4=M2=0, F6=M3=0; 
F(izl) = 0; %Set force components to zero
iecF = 2*necF-1;
F(iecF) = ecF; %Applied concentrated forces at nodes
iecM = 2*necM;
F(iecM) = ecM; %Applied concentrated moments at nodes
fprintf('Global force vector for entire system\n')
F
end

function B = apply_distributed_load(edF, dF, L)
B = sym(zeros(2*(length(L)+1), 1));
for e = 1:length(edF)
    idF = 2*e-1:2*e+2; %indices of nodes which pick up the distributed force
    B(idF) = B(idF) + 1/2*dF(e)*L(e)*[1; L(e)/6; 1; -L(e)/6];
end
end

function [Kr, Fr, Br, iR, iud] = part_system(DOF, nYr, nMr, K, F, B)
%Extract from larger system the part involving unkown generalized displacements only
iR = [2*nYr-1, 2*nMr]; %indicies of unkown reaction components (F1=V1, F5=V3) and imposed displacements (u1=v1=u5=v3=0) 
iud = setdiff(1:DOF, iR); %indices of unknown generalized displacement components
Kr = K(iud,iud) %reduced stiffness matrix
Fr = F(iud) %correspondingly
Br = B(iud) %correspondingly
end

function u = assemble_global_displacement_vector(DOF, iR, iud, ur)
u = sym('u', [1, DOF]).'; %Symbolic array of generalized coordinates: [u1, u2, ..., u6], i.e., [v1, phi1, v2, phi2, v3, phi3]
u(iR) = 0;
u(iud) = ur;
end

function F = solve_for_reactions(K, u, B, nYr, nMr)
FplusB = simplify(K*u);
F = FplusB - B;
disp('Rections:')
for i = 1:length(nYr)
    fprintf('Yr(%i) = %s\n', nYr(i), char(F(2*nYr(i)-1)))
end
for i = 1:length(nMr)
    fprintf('Mr(%i) = %s\n', nMr(i), char(F(2*nMr(i))))
end
end

function visualize_deformation(num_elements, L0, u0, undeformedObj)
%Interpolating polynomial and visualizing deformed structure
N = @(z,L) [1-3*(z/L).^2+2*(z/L).^3; z-2*z.^2/L+z.^3/L^2; 3*(z/L).^2-2*(z/L).^3; -z.^2/L+z.^3/L^2]'; %Shape functions
offset = 0;
sf = 100; %Scale factor for visualization (amplifies deformation)
for i=1:num_elements
    z = linspace(offset, offset+L0(i), 100);
    v = @(z) N(z-z(1),L0(i))*u0(2*i-1:2*i+2); %Shape functions receive local z coord for element
    v = v(z);
    deformedObj = plot(z, sf*v, 'r-', 'linewidth', 2); %positive deflection downward
    offset = offset + L0(i);
end
legend([undeformedObj, deformedObj], {'undeformed','deformed'})
% Scale down the y-tick values
yticks = get(gca, 'YTick'); % Get current y-ticks
yticklabels = yticks/sf;  % Scale down the tick values
set(gca, 'YTickLabel', yticklabels); % Scale down the y-tick values by a factor of sf
end