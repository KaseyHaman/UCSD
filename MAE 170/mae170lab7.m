% clc;
% clear all;
% close all;
% 
% load('Part2Data.mat')
% 
% clc;
% % Constants
% speed_of_sound = 343;  % Speed of sound in air (m/s)
% 
% %% 1)
% % Extract relevant data from recMatrix_sig
% time_delays = zeros(1, pointsx);  % Initialize time delay array
% x_positions = (0:pointsx-1) * move / 1000;  % Convert x positions to meters
% 
% for i = 1:pointsx
%     % Find the peak of the signal at the central y position
%     [~, idx(i)] = max(recMatrix_sig(:, i, ceil(pointsy/2)));
%     [~, idx2(i)] = max(recMatrix_ref(:, i, ceil(pointsy/2)));
%     time_delays(i) = (idx2(i)*dt)-(idx(i)*dt);
% end
% 
% for i = 1:pointsx
%     % Find the peak of the signal at the central y position
%     
%     %time_delays(i) = t(idx2) * 1e3;  % Convert time to milliseconds
% end
% 
% % Define errors
% x_error = move / 1000 * 0.05;  % Assume 5% error in position
% time_error = 0.1;  % Assume a fixed time error of 0.1 milliseconds
% 
% % Create error vectors
% x_neg = x_error * ones(size(x_positions));
% x_pos = x_error * ones(size(x_positions));
% y_neg = time_error * ones(size(time_delays));
% y_pos = time_error * ones(size(time_delays));
% 
% % Calculate theoretical time delays
% theoretical_time_delays = ((0.35-x_positions) / speed_of_sound) * 1e3;  % Convert to milliseconds
% 
% % Perform linear regression to determine experimental sound speed
% p = polyfit(x_positions, time_delays, 1);
% experimental_speed_of_sound = 1 / (p(1) / 1e3);  % Convert to m/s from ms/m
% 
% % Generate fit line
% fit_time_delays = polyval(p, x_positions);
% 
% % Plot with error bars
% figure;
% errorbar(x_positions, time_delays, y_neg, y_pos, x_neg, x_pos, 'o');
% hold on;
% plot(x_positions, theoretical_time_delays, '-r', 'LineWidth', 2);  % Theoretical line
% plot(x_positions, fit_time_delays, '--b', 'LineWidth', 2);  % Fitted line
% xlabel('Position (meters)');
% ylabel('Time Delay (milliseconds)');
% title('Time Delay of Sound Pulse Arrival along the X-Axis');
% legend('Data with error bars', 'Theoretical line (343 m/s)', 'Fitted line');
% grid on;
% set(gca, 'FontSize', 14, 'LineWidth', 2);
% 
% % Display experimental speed of sound
% disp(['Experimental speed of sound: ', num2str(experimental_speed_of_sound), ' m/s']);
% 
% %% 2)
% % Normalize the signal amplitudes
% normalized_recMatrix_sig = recMatrix_sig ./ max(abs(recMatrix_sig), [], 'all');
% 
% % Calculate time in milliseconds
% time_ms = t * 1e3;
% 
% % Create the pcolor plot
% figure;
% pcolor(x_positions, time_ms, normalized_recMatrix_sig(:, :, ceil(pointsy/2)));
% shading flat;
% xlabel('Position along X-axis (m)');
% ylabel('Time (ms)');
% title('Normalized Signal Amplitude vs. Position and Time');
% colorbar;
% c = colorbar;
% set(c, 'FontSize', 20);
% 
% % Overlay the plot with a line corresponding to the speed of sound in air
% hold on;
% sound_speed_line = (speed_of_sound * 1e3/ 0.35-x_positions) ;  % Convert to milliseconds
% plot(x_positions, sound_speed_line, '-r', 'LineWidth', 2);
% legend('Normalized Signal Amplitude', 'Speed of Sound Line');
% grid on;
% set(gca, 'FontSize', 14, 'LineWidth', 2);


