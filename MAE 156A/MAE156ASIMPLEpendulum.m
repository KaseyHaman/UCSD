clear all;
close all;
clc;
format long;
name = 'Kasey Haman';
id = 'A16978114';

%% Simple Pendulum Oscillation Without Friction



%Given Variables
g = 9.81; % m/s^2
L = 0.154; % Length to center of mass in meters
m = 0.3123; %total mass in kg
I_total = 0.00837; % Intertia in kg*m^2
theta0 = 30; % Theta initial condition
theta0_rad = deg2rad(theta0); % Theta initial condition in rad
omega0 = 0; % Theta dot initial condition
tspan = [0, 10]; % Set time span in seconds
dt = 0.001; % Set time step in seconds
radii = [1/8, 1/4, 1/2] * 0.0254; % Seperate radius values in meters

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

    % RK4 Approximation:
    for j = 1:length(t)-1
        % Calculate RK4 Constants for both theta dot and theta with derived
        % equation
        k1 = dt * [omega(j); (-m*g * L * (theta(j))) / I_total];
        k2 = dt * [omega(j) + k1(2)/2  ; (-m*g * L * (theta(j) + k1(1)/2)) / I_total];
        k3 = dt * [omega(j) + k2(2)/2  ; (-m*g * L * (theta(j) + k2(1)/2)) / I_total];
        k4 = dt * [omega(j) + k3(2)  ; (-m*g * L * (theta(j) + k3(1))) / I_total];

        % Store theta and omega
        theta(j+1) = theta(j) + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1)) / 6;
        omega(j+1) = omega(j) + (k1(2) + 2*k2(2) + 2*k3(2) + k4(2)) / 6;


        % Kinetic Energy: (1/2) * I_total * omega^2
        KE(j) = (0.5 * I_total * omega(j)^2);
        % Potential Energy: m_total * g * L * (1 - cos(theta))
        PE(j) = m * g * abs(L) * (1-cos(theta(j)));
        % Total Energy: KE + PE
        TE(j) = KE(j) + PE(j);
    end

    % Store results
    results = theta; % Store angles
    time = t; % Store time

    %Energy Results
    PEresults = PE; 
    KEresults = KE; 
    TEresults = TE; 

    %Theoretical results
    omega_n = sqrt(m*g*L/I_total); %Derivation for omega
    thetaCFS = theta0_rad * cos(omega_n * t);  %Closed form theta soln


% Plot results
figure;
    plot(time, rad2deg(results), 'r', 'LineWidth', 1.5);
    hold on;
    plot(time, rad2deg(thetaCFS), 'b--', 'LineWidth', 3)
    xlabel('Time (s)');
    ylabel('\theta (degrees)');
    title(['Pendulum Oscillation']);
    axis([0 10 -40 40])
    legend('ODE Simulation', 'Closed Form Solution', 'Location', 'best');
    grid on;



% Plot energy (KE, PE, TE)
figure;
    plot(time, KEresults, 'r', 'LineWidth', 1.5); hold on;
    plot(time, PEresults, 'b', 'LineWidth', 1.5);
    plot(time, TEresults, 'k', 'LineWidth', 0.5);
    xlabel('Time (s)');
    ylabel('Energy (J)');
    title(['Energy']);
    axis([0 10 -0.015 0.08])
    legend('KE (Red)', 'PE (Blue)', 'TE (Black)', 'Location', 'best');
    grid on;

