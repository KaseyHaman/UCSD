clc;clear all;

%% Part 1
clc;clear all;
%Load data
load('cylinder_Re100.mat')


% udata = u(1, :, :);
% udata = squeeze(udata);

%Size u
[a b c] = size(u);

%Run simulation
for t = 1:a
figure(1);
    ustep = squeeze(u(t, :, :));
    % Subplot 1
    subplot(1, 2, 1);
    pcolor(x, y, ustep);
    shading interp;
    axis equal tight;
    title('u');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    drawnow;

    vstep = squeeze(v(t, :, :));
    % Subplot 1
    subplot(1, 2, 2);
    pcolor(x, y, vstep);
    shading interp;
    axis equal tight;
    title('v');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    drawnow;
end


%% Part 2
clc;clear all;
%Load data
load('cylinder_Re100.mat')
%Size u
[a b c] = size(u);

u2 = squeeze(mean(u(150:end, :, :), 1));
v2 = squeeze(mean(v(150:end, :, :), 1));
figure(2);

  % Subplot 1
    subplot(1, 2, 1);
    pcolor(x, y, u2);
    shading interp;
    axis equal tight;
    title('u');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    % Subplot 2
    subplot(1, 2, 2);
    pcolor(x, y, v2);
    shading interp;
    axis equal tight;
    title('v');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);


%% Part 3
clc;clear all;
%Load data
load('cylinder_Re100.mat')

%Plot part 2

%Size u
[a b c] = size(u);

u2 = squeeze(mean(u(150:end, :, :), 1));
v2 = squeeze(mean(v(150:end, :, :), 1));

figure(3);

    pcolor(x, y, u2);
    shading interp;
    axis equal tight;
    title('u');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);


    hold on;

%Run simulation for streamline
streamline(x', y', u2', v2', -4*ones(10,1),[-4 -3 -2 -1 -0.01 0.01 1 2 3 4]);


%% Part 4

clc;clear all;
%Load data
load('cylinder_Re100.mat')


% udata = u(1, :, :);
% udata = squeeze(udata);

%Size u
[a b c] = size(u);

u2 = squeeze(mean(u(150:end, :, :), 1));
v2 = squeeze(mean(v(150:end, :, :), 1));


%Run simulation
for t = 1:a
figure(1);
    ustep = squeeze(u(t, :, :));
    uprime = ustep - u2;

    % Subplot 1
    subplot(1, 2, 1);
    pcolor(x, y, uprime);
    shading interp;
    axis equal tight;
    title('u');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    drawnow;

    vstep = squeeze(v(t, :, :));
    vprime = vstep - v2;

    % Subplot 1
    subplot(1, 2, 2);
    pcolor(x, y, vprime);
    shading interp;
    axis equal tight;
    title('v');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    drawnow;
end


%% Part 5
clc;clear all;
%Load data
load('cylinder_Re100.mat')


% udata = u(1, :, :);
% udata = squeeze(udata);

%Size u
[a b c] = size(u);

u2 = squeeze(mean(u(150:end, :, :), 1));
v2 = squeeze(mean(v(150:end, :, :), 1));

uprime = zeros(150, b, c);
vprime = zeros(150, b, c);

for t = 150:a

   ustep = squeeze(u(t, :, :));
   uprime(t, :, :) = ustep - u2;

    vstep = squeeze(v(t, :, :));
    vprime(t, :, :) = vstep - v2;

end

figure(1);

    k = 1/2*(mean(uprime.^2)+mean(vprime.^2));

    %Plot
    pcolor(x, y, squeeze(k));
    shading interp;
    axis equal tight;
    title('u');
    xlabel('x');
    ylabel('y');
    rectangle('Position', [-0.5 -0.5 1 1], 'Curvature', [1 1], ...
              'LineStyle', 'none', 'FaceColor', [1 1 1]);
    colorbar
    drawnow;
