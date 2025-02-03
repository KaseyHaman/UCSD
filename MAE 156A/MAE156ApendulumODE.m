clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%%MAE 156A Pendulum ODE Problem
g = 9.81; % m/s^2
L = 0.154; % Length to center of mass in meters
m = 0.3123; %total mass in kg
I_total = 0.00839; % Intertia in kg*m^2
m_acrylic = 0.0923; % Mass of acrylic (kg)
m_bolt = 0.055;    % Mass of each bolt (kg)
Lzz = 4.88*10^-4;  %Inertia of acrylic about COM in kg*m^2
Lacrylic = 0.1152; %Length from center of rotation to acrylic COM in m
theta0 = 30; % Theta initial condition
theta0_rad = deg2rad(theta0); % Theta initial condition in rad
omega0 = 0; % Theta dot initial condition
tspan = [0, 100]; % Set time span in seconds
dt = 0.001; % Set time step in seconds
radii = [1/8, 1/4, 1/2] * 1/2 * 0.0254; % Seperate radius values in meters
mu = 0.11; % Coefficient of friction


%Calculate Inertia
bolt_coords = [
    0, 5.45 * 0.0254;        % Bolt 1 location in meters
    0, (5.45 + 2.58) * 0.0254; % Bolt 2 location in meters
    2.58 / 2 * 0.0254, (5.45 + 2.58 / 2) * 0.0254; % Bolt 3 location in meters
    -2.58 / 2 * 0.0254, (5.45 + 2.58 / 2) * 0.0254  % Bolt 4 location in meters
];

% Moment of inertia of acrylic
I_acrylic = Lzz + m_acrylic * Lacrylic^2; %Using given values for acrylic in metric

% Calculate moment of inertia for each bolt using parallel axis theorem
I_bolts = 0;
for i = 1:size(bolt_coords, 1)
    r = sqrt(bolt_coords(i, 1)^2 + bolt_coords(i, 2)^2); % distance to z-axis
    I_bolts = I_bolts + m_bolt * r^2;
end

% Total moment of inertia (sum of acrylic and bolts)
I_total = I_acrylic + I_bolts;

%Total mass
m = 4*m_bolt + m_acrylic;

% Loop through each radius (bushing size)
for i = 1:length(radii)
    r = radii(i); % Current bushing radius

    y0 = [theta0_rad; omega0]; % Set IC's
    t = tspan(1):dt:tspan(2); % Set time vector

    % Pre-allocate for storing results
    theta = zeros(size(t)); 
    omega = zeros(size(t)); 
    theta(1) = y0(1); 
    omega(1) = y0(2); 
    KE = zeros(size(t)); 
    PE = zeros(size(t)); 
    TE = zeros(size(t)); 


    % Perform RK4 integration for each time step
    for j = 1:length(t)-1
        % Calculate the four slopes (k1, k2, k3, k4) for RK4
        k1 = dt * pendulumODE([theta(j); omega(j)], r, g, L, m, I_total, mu);
        k2 = dt * pendulumODE([theta(j) + k1(1)/2; omega(j) + k1(2)/2], r, g, L, m, I_total, mu);
        k3 = dt * pendulumODE([theta(j) + k2(1)/2; omega(j) + k2(2)/2], r, g, L, m, I_total, mu);
        k4 = dt * pendulumODE([theta(j) + k3(1); omega(j) + k3(2)], r, g, L, m, I_total, mu);

        % Update the state (theta, omega)
        theta(j+1) = theta(j) + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1)) / 6;
        omega(j+1) = omega(j) + (k1(2) + 2*k2(2) + 2*k3(2) + k4(2)) / 6;


        % Kinetic Energy: (1/2) * I_total * omega^2
        KE(j) = (0.5 * I_total * omega(j)^2);
        % Potential Energy: m_total * g * L * (1 - cos(theta))
        PE(j) = m * g * abs(L) * (1-cos(theta(j)));
        % Total Energy: KE + PE
        TE(j) = KE(j) + PE(j);

   %Find and mark the end of the oscillation
   if abs(rad2deg(theta(j+1))) < 0.1 && abs(omega(j+1)) < 0.1
    Diameter_inch = r*2 / 0.0254; 
    disp(['Pendulum with diameter ', num2str(Diameter_inch), ' inch stopped oscillating at t = ', num2str(t(j)), ' seconds.']);
        stoptime(i) = t(j);
        totalt{i} = t;
        stoptheta(i) = theta(j);
        break; 
   end

    end

    % Store results
    results{i} = theta; % Store angles
    time{i} = t; % Store time

    %Energy Results
    PEresults{i} = PE; 
    KEresults{i} = KE; 
    TEresults{i} = TE; 

end

% Plot results in subplots
figure;
for i = 1:length(radii)
    subplot(3, 1, i);
    plot(time{i}, rad2deg(results{i}), 'LineWidth', 1.5);
    hold on;
    plot(stoptime(i), stoptheta(i), 'o', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (degrees)');
    title(['Pendulum Oscillation (diameter = ', num2str(radii(i) / (1/2 * 0.0254)), ' inch)']);
    legend('Oscillation', ['Stopped at t = ', num2str(stoptime(i), '%.2f'), ' s'], 'Location', 'northeast');
    grid on;
end


% Plot energy (KE, PE, TE) for all radii
figure;
for i = 1:length(radii)
    subplot(3, 1, i);
    plot(time{i}, KEresults{i}, 'r', 'LineWidth', 1.5); hold on;
    plot(time{i}, PEresults{i}, 'b', 'LineWidth', 1.5);
    plot(time{i}, TEresults{i}, 'k', 'LineWidth', 0.5);
    xline(stoptime(i), '--k', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Energy (J)');
    title(['Energy (diameter = ', num2str(radii(i) / (1/2 * 0.0254)), ' inch)']);
    legend('KE (Red)', 'PE (Blue)', 'TE (Black)', 'Stop of Oscillation', 'Location', 'southeast');
    grid on;
end


% Function to calculate the ODE 
function pendulumSOLN = pendulumODE(y, r, g, L, m, I, mu)

    % Angular acceleration (theta double dot)
    pendulumSOLN = [
        y(2); % theta'
        -m*g*sin(y(1))*L/I + calculateFriction(y(1), y(2), mu, m, g, r, I, L)
    ];
end

% Function to calculate friction force
function friction = calculateFriction(theta, omega, mu, m, g, r, I, L)
    if omega > 0
        friction = (-r*m*mu*(g*cos(theta) + omega^2*L)) / I; % Opposes motion when omega > 0
    elseif omega < 0
        friction = (r*m*mu*(g*cos(theta) + omega^2*L)) / I; % Opposes motion when omega < 0
    else
        friction = 0; % No friction when omega = 0
    end
end




