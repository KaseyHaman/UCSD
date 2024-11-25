%% Part 1 Code
load("acousticscan739394385.mat")

%Question 1
%Note: recMatrix_sig contains running average of each signal as it is
%calculated and wave_sig is then the last signal measured and used in the
%running average

figure(1)
hold on
title("Last Acquired Signal vs Averaged Signal")
plot(t,wave_sig)
plot(t,recMatrix_sig)
xlabel("Time (s)")
ylabel("Amplitude (V)")
legend("Last Acquired Signal", "Averaged Signal")

%Question 2
figure(2)
hold on
plot(t, recMatrix_sig)
plot(t, recMatrix_ref)

%% Part 2 Code


% load("Part2Data.mat")
% timedelay = zeros(30,15);
% for i = 1:30
%    
%         findpeaks(recMatrix_sig(:,i,1),MinPeakDistance=12);
%         [pk2, lc2] = findpeaks(recMatrix_ref(:,i,1),MinPeakProminence=1);
%         actualpulse = lc2(1);
%         measuredpulse = input('Input the position of first pulse');
%         timedelay(i,1) = t(measuredpulse)-t(actualpulse);
% 
% 
% end

% xpos = zeros(1,30) +0.01;
% xneg = zeros(1,30) +0.01;
% ypos = zeros(1,30) + 0.000125;
% yneg = zeros(1,30) + 0.000125;
% delay = timedelay(:,1);
% figure(1)
% hold on, box on, grid on
% 
% errorbar(x/1000,timedelay(:,1),ypos,yneg,xpos,xneg)
% p = polyfit(x/1000,timedelay(:,1),1);
% fit = polyval(p,x/1000);
% plot(x/1000,fit)

%% Alternative plot
xpos = zeros(1,30) +0.001;
xneg = zeros(1,30) +0.001;
ypos = zeros(1,30) + 0.000125;
yneg = zeros(1,30) + 0.000125;
delay = timedelay(:,1);
figure(2)
hold on,box on, grid on
distance = sqrt((0.35-(x/1000)).^2 + 0.15^2);
errorbar(distance,delay,ypos,yneg,xpos,xneg,'o')
p = polyfit(distance,delay,1);
fit = polyval(p,distance);
plot(distance,fit,'LineWidth',2)
title('Time Delay vs Distance from Microphone')
xlabel('Distance (m)')
ylabel('Time Delay (s)')
y = 1/340*distance + p(2);
plot(distance,y)
legend('Errorbar', 'Experimental Speed of sound','Actual Speed of Sound')


speedsound = 1/p(1)