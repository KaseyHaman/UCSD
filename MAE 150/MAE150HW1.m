clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1%

%Part a
star = [4, 0;
        1, 1;
        0, 7;
       -1, 1;
       -4, 0;
       -1, -1;
        0, -7;
        1, -1;
        4, 0]; 

plot(star(:,1), star(:,2), '-b', 'LineWidth', 1)
hold on;
axis([-30 30 -30 30]);
title('Star Shaped Transformations');
grid on;
axis equal;

%Part b
plot([-40 40], [0, 0], 'b-', 'LineWidth', 0.5)
plot([0 0], [-40, 40], 'b-', 'LineWidth', 0.5)
xlabel('x');
ylabel('y');


%Part c
theta = pi/3;
rot_matrix = [cos(theta), -sin(theta), 0;
    sin(theta), cos(theta) 0;
    0, 0, 1];

%Modify Star Matrix
augmented_star = [star, ones(size(star, 1), 1)]';

rot_star = rot_matrix * augmented_star;
plot(rot_star(1, :), rot_star(2,:), 'r-', 'LineWidth', 1)

%Part d
translation_matrix = [1 0 7; 0 1 3; 0 0 1];
translated_star = translation_matrix * rot_star;
plot(translated_star(1, :), translated_star(2,:), 'g-', 'LineWidth', 1)

%Part e
scaling_matrix = [2 0 0; 0 4 0; 0 0 1];
scaled_star = scaling_matrix * translated_star;
plot(scaled_star(1, :), scaled_star(2,:), 'm-', 'LineWidth', 1);

%Part f
reflection_matrix = [1 0 0; 0 -1 0; 0 0 1];
reflected_star = reflection_matrix * scaled_star;
plot(reflected_star(1, :), reflected_star(2,:), 'c-', 'LineWidth', 1)


%%Problem 2
shape = [3 6; 6 6; 6 7; 4 7; 4 11; 3 11];
plot([shape(:,1); shape(1,1)], [shape(:,2); shape(1,2)], 'k-', 'LineWidth', 2);
LinePoint1 = [13, 5];
LinePoint2 = [-35, 2];

hold on;
title('Reflection Across Arbitrary Line')
plot([13, -35], [5, 2], '-b', 'LineWidth', 1)
plot(LinePoint1(1), LinePoint1(2), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 5); % Point A
plot(LinePoint2(1), LinePoint2(2), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 5); % Point B
axis equal;

%Translation Matrix
translation_matrix = [1 0 -LinePoint1(1); 0 1 -LinePoint1(2); 0 0 1];
invtranslation_matrix = [1 0 LinePoint1(1); 0 1 LinePoint1(2); 0 0 1];


%Find angle of arbitrary line
theta = atan((LinePoint2(2) - LinePoint1(2))/(LinePoint2(1) - LinePoint1(1)));

%Rotate line to axis (CLOCKWISE ROTATION)
rot_matrix = [cos(-theta), -sin(-theta), 0;
    sin(-theta), cos(-theta) 0;
    0, 0, 1];

%Reflect across axis
reflection_matrix = [1 0 0; 0 -1 0; 0 0 1];

%Rotate back the axis (COUNTER-CLOCKWISE ROTATION)
invrot_matrix = [cos(theta), -sin(theta), 0;
    sin(theta), cos(theta) 0;
    0, 0, 1];

%Prepare shape for multiplication
augmented_shape = [shape'; ones(1, size(shape, 1))];

%Multiply matricies in proper order
reflected_shape = (invtranslation_matrix * invrot_matrix * reflection_matrix * rot_matrix * translation_matrix * augmented_shape)';

%Plot reflected shape
plot([reflected_shape(:,1); reflected_shape(1,1)], [reflected_shape(:,2); reflected_shape(1,2)], 'k--', 'LineWidth', 2);


%%Problem 3
%Part a
rectangle = [0, 0; 4, 0; 4, 1; 0,1; 0 0];
plot(rectangle(:, 1), rectangle(:, 2), '-b', 'LineWidth', 1)
hold on;
axis([-25 25 -25 25])
xlabel('x');
ylabel('y');
axis equal;
grid on;
title('Rotation Animation')

%Part b
B = [10 0];
plot(B(1), B(2), 'ro', 'LineWidth', 1)

%Part c
theta = 6*pi/100;

translation_matrix = [1 0 -B(1); 0 1 -B(2); 0 0 1];

invtranslation_matrix = [1 0 B(1); 0 1 B(2); 0 0 1];

rot_matrix = [cos(theta), -sin(theta), 0;
    sin(theta), cos(theta) 0;
    0, 0, 1];

augmented_rectangle = [rectangle, ones(size(rectangle, 1) ,1)]';

rot_rectangle = invtranslation_matrix * rot_matrix * translation_matrix * augmented_rectangle;

%Part d
for n=1:100

    clf;

    theta = 6*pi*n/100;

rot_matrix = [cos(theta), -sin(theta), 0;
    sin(theta), cos(theta) 0;
    0, 0, 1];

rot_rectangle = invtranslation_matrix * rot_matrix * translation_matrix * augmented_rectangle;

plot(rot_rectangle(1, :), rot_rectangle(2, :), '-b', 'LineWidth', 1)
axis([-25 25 -25 25])
xlabel('x');
ylabel('y');
axis equal;
grid on;
title('Rotation Animation')
hold on;
plot(B(1), B(2), 'ro', 'LineWidth', 1)

pause(0.005)
end


%For future reference:

%Rotation Matrix (CCW)
rot_matrix = [cos(theta), -sin(theta), 0;
    sin(theta), cos(theta) 0;
    0, 0, 1];
%Rotate (CLOCKWISE ROTATION)
rot_matrix = [cos(-theta), -sin(-theta), 0;
    sin(-theta), cos(-theta) 0;
    0, 0, 1];

%Translation Matrix
translation_matrix = [1 0 -LinePoint1(1); 0 1 -LinePoint1(2); 0 0 1];
invtranslation_matrix = [1 0 LinePoint1(1); 0 1 LinePoint1(2); 0 0 1];

%Scaling Matrix
scaling_matrix = [2 0 0; 0 4 0; 0 0 1];


%Reflection Matrix
reflection_matrix = [1 0 0; 0 -1 0; 0 0 1];





