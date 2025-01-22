clear all; clc; close all;

%Truss member specifications
E = 200e9; % Young's modulus (Pa)
A = 0.001; % cross-sectional area (m^2)
L = 1; % length (m)

%Angle each element makes with horizontal (measured in the CCW direction)
phi = [0,  0, 60, 180, 180, 60, 120, 60, 120, 60, 120]*pi/180;

%Specify nodes of each element
iL = [1,  2,  3,  4,  5,  6,  5,  1,  2,  2,  3]; %nodes having lower index number of the tree elements
iH = [2,  3,  4,  5,  6,  7,  7,  6,  6,  5,  5]; %nodes having higher index number of the tree elements

%z and y coordinates of nodes 1, 2, etc.
z = L*[-2, -1,   0,   1/2,       -1/2,      -3/2,      -1];
y = L*[ 0,  0,   0, sqrt(3)/2, sqrt(3)/2, sqrt(3)/2, sqrt(3)];

%List of node coordinates for visualizing all the bar elements
Z = [z(1), z(2), z(3), z(4), z(5), z(6), z(1), z(2), z(5), z(3), z(2), z(6), z(7), z(5)];
Y = [y(1), y(2), y(3), y(4), y(5), y(6), y(1), y(2), y(5), y(3), y(2), y(6), y(7), y(5)];

% Visualize the undeformed structure
figure(1)
set(0,'defaultAxesFontName', 'Times')
undeformedObj = line(Z,Y,'color','k','linewidth',2);
axis([-2 0.5 -0.5 2])
axis equal
grid on
title('Undeformed and deformed truss') 
set(gca,'fontsize',14);
xlabel('{\itz} (m)');
ylabel('{\ity} (m)');
hold on

%Visualize supports
plot(0,-0.07,'b^','MarkerSize',10)
plot([-0.1 0.1],[-0.11 -0.11],'b-')
plot(0.55,0.87,'bo','MarkerSize',10)
plot([0.6 0.6],[0.80 0.94],'b-')

%Visualize loads
quiver(-2, 0, 0, -0.4, 'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 5); % Disable scaling
quiver(-1, 1.74, -0.4, 0, 'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 5); % Disable scaling

pause(0.05)

s = 100; %Scale factor

%Definition of the element stiffness matrices
num_nodes = length(z);
num_elements = length(iL);
DOF = 2*num_nodes; %number of degrees of freedom (2 per node)
K = zeros(DOF);
for i=1:num_elements
    Keg = zeros(DOF); %Preallocation of element stiffness matrix expanded to global size
    sprintf('Individual stiffness matrix for element %d', i)
    s2 = sin(phi(i))^2;
    c2 = cos(phi(i))^2;
    sc = sin(phi(i))*cos(phi(i));
    Ke = E*A/L*[c2,  sc, -c2, -sc;
                   sc,  s2, -sc, -s2;
                  -c2, -sc,  c2,  sc;
                  -sc, -s2,  sc,  s2];
    Keg([2*iL(i)-1, 2*iL(i), 2*iH(i)-1, 2*iH(i)],[2*iL(i)-1, 2*iL(i), 2*iH(i)-1, 2*iH(i)]) = Ke;
    K = K + Keg;
end

%Apply external loading conditions in z and y direction per node
F = zeros(DOF,1);
F(2) = -10e3; %kN (=F1y)
F(DOF-1) = -10e3; %kN (=F7z)

%Extract from larger system the part involving unkown displacements only
inds = setdiff(1:DOF, 5:7); %since u5=z3=0, u6=y3=0, u7=z4=0 
Kr = K(inds,inds);
Fr = F(inds); 

%Invert reduced matrix and solve for unknown displacements
ur = inv(Kr)*Fr;

%Construct full global displacement vector
u = zeros(DOF,1);
u(inds) = ur

%Solve for the unknown reactions
F = K*u

%Calculate coordinates of deformed structures while applying magnification factor s 
zd = z + s*u([1:2:DOF])';
yd = y + s*u([2:2:DOF])';

%List of node coordinates for visualizing all the deformed bar elements
Zd = [zd(1), zd(2), zd(3), zd(4), zd(5), zd(6), zd(1), zd(2), zd(5), zd(3), zd(2), zd(6), zd(7), zd(5)];
Yd = [yd(1), yd(2), yd(3), yd(4), yd(5), yd(6), yd(1), yd(2), yd(5), yd(3), yd(2), yd(6), yd(7), yd(5)];

%Visualize deformed sructure
deformedObj = line(Zd,Yd,'color','r','linewidth',2);
legend([undeformedObj, deformedObj],{'undeformed','deformed'})
