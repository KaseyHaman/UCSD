clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1
%Part A
theta = pi/6;
shear_matrix = [1 tan(theta) 0; 0 1 0; 0 0 1];

%Part B
Rectangle = [0 0;
        0, 4;
        3, 4;
       3, 0;
       0, 0];
L = [0 0; 0 5; 1 5; 1 1; 3 1; 3 0];


figure(1)
plot(Rectangle(:,1), Rectangle(:,2), '-k', 'LineWidth', 1)
hold on;
axis([0 6 0 4]);
title('Rectangle Shear');
xlabel('x');
ylabel('y');
grid on;

Aug_Rectangle = [Rectangle ones(size(Rectangle, 1), 1)]';

Shear_rectangle = (shear_matrix * Aug_Rectangle)';

plot(Shear_rectangle(:,1), Shear_rectangle(:,2), '-b', 'LineWidth', 1)


figure(2)
plot([L(:,1); L(1,1)], [L(:,2); L(1,2)], 'k-', 'LineWidth', 2);
hold on;
axis([0 6 0 5]);
title('L Shape Shear');
xlabel('x');
ylabel('y');
grid on;


Aug_L = [L ones(size(L, 1), 1)]';
Shear_L = (shear_matrix * Aug_L)';

plot(Shear_L(:,1), Shear_L(:,2), '-b', 'LineWidth', 1)

%Part C

invshear_matrix = [1 -tan(theta) 0; 0 1 0; 0 0 1];

A = invshear_matrix * shear_matrix; %Proof

%Part D
%Based on our notes, a matrix is orthogonal if S^T*S = I
B = shear_matrix' * shear_matrix;
disp(B)
%Because B does not equal the identity matrix, it is not orthogonal.

%Part E
DetShearingMatrix = det(shear_matrix);
disp(DetShearingMatrix)
%As stated in the notes, the determinant representes the volume
%(or area) scaling factor induced by the given mapping. So, because the det = 1,
% the area of the figure is preserved in the transformation.


