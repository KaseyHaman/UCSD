clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%Problem 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part a
% Coordinates of the star-shaped object
star_points = [4, 0;
               1, 1;
               0, 7;
               -1, 1;
               -4, 0;
               -1, -1;
               0, -7;
               1, -1;
               4, 0]; % Close the shape

% Plot the star
figure;
plot(star_points(:,1), star_points(:,2), 'b-', 'LineWidth', 2);
hold on;

% Set the axis limits and labels
axis([-30 30 -30 30]);
xlabel('x');
ylabel('y');
axis equal;
grid on;
title('Star-shaped object and transformations');


%Part b
% Plot the x and y axes
plot([-30, 30], [0, 0], 'c-', 'LineWidth', 1); % x-axis
plot([0, 0], [-30, 30], 'c-', 'LineWidth', 1); % y-axis

%Part c
% Define the rotation angle (60 degrees in radians)
theta = deg2rad(60);

% Homogeneous rotation matrix
rotation_matrix = [cos(theta), -sin(theta), 0;
                   sin(theta), cos(theta), 0;
                   0, 0, 1];

% Display the rotation matrix
disp('Homogeneous rotation matrix (60 degrees):');
disp(rotation_matrix);

% Apply the rotation matrix to the star points
homogeneous_star_points = [star_points, ones(size(star_points, 1), 1)]'; % Convert to homogeneous coordinates
rotated_star_points = (rotation_matrix * homogeneous_star_points)'; % Apply rotation

% Plot the rotated star using red lines
plot(rotated_star_points(:,1), rotated_star_points(:,2), 'r-', 'LineWidth', 2);

%Part d
% Translation vector
t_x = 7;
t_y = 3;

% Homogeneous translation matrix
translation_matrix = [1, 0, t_x;
                      0, 1, t_y;
                      0, 0, 1];

% Display the translation matrix
disp('Homogeneous translation matrix:');
disp(translation_matrix);

% Apply the translation matrix to the star points
translated_star_points = (translation_matrix * rotated_star_points')'; % Apply translation

% Plot the translated star using green lines
plot(translated_star_points(:,1), translated_star_points(:,2), 'g-', 'LineWidth', 2);


%Part e
% Scaling factors
s_x = 2;
s_y = 4;

% Homogeneous scaling matrix
scaling_matrix = [s_x, 0, 0;
                  0, s_y, 0;
                  0, 0, 1];

% Display the scaling matrix
disp('Homogeneous scaling matrix:');
disp(scaling_matrix);

% Apply the scaling matrix to the star points
scaled_star_points = (scaling_matrix * translated_star_points')'; % Apply scaling

% Plot the scaled star using magenta lines
plot(scaled_star_points(:,1), scaled_star_points(:,2), 'm-', 'LineWidth', 2);


%Part f
% Homogeneous reflection matrix (over x-axis)
reflection_matrix_x = [1, 0, 0;
                       0, -1, 0;
                       0, 0, 1];

% Display the reflection matrix
disp('Homogeneous reflection matrix (over x-axis):');
disp(reflection_matrix_x);

% Apply the reflection matrix to the star points
reflected_star_points = (reflection_matrix_x * scaled_star_points')'; % Apply reflection

% Plot the reflected star using cyan lines
plot(reflected_star_points(:,1), reflected_star_points(:,2), 'c-', 'LineWidth', 2);



%%Problem 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the coordinates of the object (L-shaped object as an example)
object_points = [3, 6;
                 6, 6;
                 6, 7;
                 4, 7;
                 4, 11;
                 3, 11];

% Define the reflection line by two points, A and B
A = [13, 5];
B = [-35, 2];

% Plot the original object using solid lines
figure;
plot([object_points(:,1); object_points(1,1)], [object_points(:,2); object_points(1,2)], 'b-', 'LineWidth', 2);
hold on;

% Plot the reflection line
plot([A(1), B(1)], [A(2), B(2)], 'k--', 'LineWidth', 1.5);
plot(A(1), A(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Point A
plot(B(1), B(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Point B

% Set axis limits and labels
axis equal;
xlabel('x');
ylabel('y');
grid on;
title('Reflection of Object Across Line');
hold on;

% Homogeneous coordinates for the object points
homogeneous_object_points = [object_points, ones(size(object_points, 1), 1)];

%% Step 1: Translation matrix (to move A to the origin)
translation_to_origin = [1, 0, -A(1);
                         0, 1, -A(2);
                         0, 0, 1];

%% Step 2: Calculate the angle of rotation to align line AB with the x-axis
delta_y = B(2) - A(2);
delta_x = B(1) - A(1);
angle = atan2(delta_y, delta_x); % Angle to rotate line to the x-axis

%% Step 3: Rotation matrix (to align line AB with the x-axis)
rotation_matrix = [cos(-angle), -sin(-angle), 0;
                   sin(-angle), cos(-angle), 0;
                   0, 0, 1];

%% Step 4: Reflection matrix (reflect across the x-axis)
reflection_matrix = [1, 0, 0;
                     0, -1, 0;
                     0, 0, 1];

%% Step 5: Inverse rotation matrix (to undo the rotation after reflection)
inverse_rotation_matrix = [cos(angle), -sin(angle), 0;
                           sin(angle), cos(angle), 0;
                           0, 0, 1];

%% Step 6: Inverse translation matrix (to undo the translation after reflection)
translation_back = [1, 0, A(1);
                    0, 1, A(2);
                    0, 0, 1];

%% Final transformation matrix
transformation_matrix = translation_back * inverse_rotation_matrix * reflection_matrix * rotation_matrix * translation_to_origin;

%% Apply the transformation matrix to the object points
reflected_object_points = (transformation_matrix * homogeneous_object_points')';

%% Plot the reflected object using dashed lines
plot([reflected_object_points(:,1); reflected_object_points(1,1)], [reflected_object_points(:,2); reflected_object_points(1,2)], 'b--', 'LineWidth', 2);

% Display the reflection line and points
legend('Original Object', 'Reflection Line', 'Point A', 'Point B', 'Reflected Object');
hold off;





%%Problem 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part a
% Coordinates of the rectangle vertices
rectangle_points = [0, 0;
                    4, 0;
                    4, 1;
                    0, 1;
                    0, 0]; % Close the rectangle

% Plot the rectangle using blue lines
figure;
plot(rectangle_points(:,1), rectangle_points(:,2), 'b-', 'LineWidth', 2);
hold on;

% Set axis limits and labels
axis([-25 25 -25 25]);
xlabel('x');
ylabel('y');
axis equal;
grid on;
title('Rectangle and rotation about point B (10, 0)');

%Part b
% Define point B
B = [10, 0];

% Plot point B as a small red circle
plot(B(1), B(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

%Part c
n = 100;
theta = 6 * pi / n; % Angle of rotation for each step

% Rotation matrix around the origin
rotation_matrix_origin = [cos(theta), -sin(theta), 0;
                          sin(theta), cos(theta), 0;
                          0, 0, 1];

% Translation matrix to move the pivot to the origin
translation_to_origin = [1, 0, -B(1);
                         0, 1, -B(2);
                         0, 0, 1];

% Translation matrix to move back to the original position
translation_back = [1, 0, B(1);
                    0, 1, B(2);
                    0, 0, 1];

% Combined rotation matrix around point B
rotation_matrix_B = translation_back * rotation_matrix_origin * translation_to_origin;

% Display the rotation matrix
disp('Homogeneous rotation matrix around point B (10, 0):');
disp(rotation_matrix_B);

%Part d
% Convert rectangle points to homogeneous coordinates (add a column of ones)
homogeneous_rectangle_points = [rectangle_points, ones(size(rectangle_points, 1), 1)];

% Animate the rotation
for i = 1:n
    % Apply the rotation matrix
    rotated_rectangle_points = (rotation_matrix_B * homogeneous_rectangle_points')';
    
    % Clear previous plot
    clf;
    
    % Plot the rotated rectangle
    plot(rotated_rectangle_points(:,1), rotated_rectangle_points(:,2), 'b-', 'LineWidth', 2);
    hold on;
    
    % Plot point B again
    plot(B(1), B(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    
    % Set axis limits and labels
    axis([-25 25 -25 25]);
    xlabel('x');
    ylabel('y');
    axis equal;
    grid on;
    
    % Pause to create animation effect
    pause(0.05);
    
    % Update the rectangle points for the next iteration
    homogeneous_rectangle_points = rotated_rectangle_points;
end

hold off;
