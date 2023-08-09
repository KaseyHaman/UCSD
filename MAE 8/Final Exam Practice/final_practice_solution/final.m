clear all; close all; clc;
name = 'AAAA AAAA';
id = 'A00000000';
hw_num = 'final';
form = 'A';
format long; format compact;

%% Problem 1: (a) 12 pts, (b,c) 4 pts, total 20 pts
theta = -10:0.1:350;
x = cosd(4*theta).*cosd(theta);
y = cosd(4*theta).*sind(theta);
z = cosd(8*theta);

counter = 0;
for n = 2:length(z)-1
    if z(n) > z(n+1) && z(n) > z(n-1)
        counter = counter  + 1;
        x_tip(counter) = x(n);
        y_tip(counter) = y(n);
        z_tip(counter) = z(n);
    end
end

figure(1); hold on; 
plot3(x,y,z,'-m','LineWidth',1); 
plot3(x_tip,y_tip,z_tip,'go','MarkerFaceColor','g','MarkerSize',10);
axis tight;
title('Three-dimensional Floral Pattern'); 
legend('Floral pattern','tip of petals');
xlabel('x'); ylabel('y'); zlabel('z');
box on; grid on; view(45,60);
set(gca,'FontSize',16);
sp1a = 'See figure 1';

dxdtheta = diff(x)./diff(theta);
sp1b = dxdtheta(theta == 90);
sp1c = sum(sqrt(diff(x).^2+diff(y).^2+diff(z).^2));

%% Problem 2: (a,b,d,e) 3 pts each, (c,f) 4 pts, total 20 pts 
rng(int8(form));
for n = 1:3 
    for m = 1:3
        cellA{n,m} = randi(100,n,m);
        structA(n).field1 = char(80+[1:(3*n)]);
        structA(n).field2 = randi(100,n,m);
    end
end
sp2a = cellA{3,2}(2,2);
sp2b = mean(mean(cellA{3,3}));
for n = 1:numel(cellA)
    max_element(n) = max(cellA{n}(:));
end
sp2c = max(max_element);

sp2d = structA(2).field1(end);
sp2e = structA(3).field2(2,2);
sp2f = structA;
for n = 1:length(sp2f)
    sp2f(n).field1(3) = lower(sp2f(n).field1(3));
end


%% Problem 3: Problem 3: (a,b) 1 pts each, (c) 3 pts, (d-i) 2 pts each, (j) 4 pts, total 20 pts
x = -20:0.1:20;
f = cos(2*x).*tanh(x/10);
sp3a = x; 
sp3b = f;
sp3c = (x(2)-x(1))*( 0.5*(f(1)+f(end)) + sum(f(2:end-1)) );

counter = 0;
for i=2:length(x)-1 
    if f(i)>f(i+1) && f(i) > f(i-1) 
        counter = counter+1;
        sp3e(counter) = x(i);
        sp3f(counter) = f(i);
    end 
end 
sp3d = counter;

counter = 0;
for i=2:length(x)-1 
    if (f(i)<0 && f(i+1) > 0) || (f(i)>=0 && f(i+1)< 0) || (f(i) == 0)     
        counter = counter+1;
        sp3h(counter) = x(i);
        sp3i(counter) = f(i);
    end 
end 
sp3g = counter;

figure(2);
plot(sp3a,sp3b,'-k',sp3e,sp3f,'rs',sp3h,sp3i,'bo');
title('Function f(x), its local maxima and zero crossing');
xlabel('x'); ylabel('f(x)'); box on; grid on;
legend('f(x)','maxima','zero crossing','location','best');
set(gca,'FontSize',14);
sp3j = 'See figure 1';


%% Problem 4: (a-c) 5 pts each, total 15 points
sp4a = 0;
for m = 1:100
    for n = 1:100
        sp4a = sp4a + 1/2^(m*n);
    end
end


day = 1;
saving = 1;
balance = 1;
while balance < 10000
    day = day + 1;
    saving = 3*saving;
    balance = balance + saving;
    %fprintf('Day: %d - saving: %d - balance: %d\n', day,saving, balance)
end
sp4b = day;

k = 0;
n = 0;
series = 4*(-1)^n/(2*n+1);
LHS = abs(pi - series);
while LHS > 1e-3
    k = k + 1;
    series = 0;
    for n = 0:k
        series = series + 4*(-1)^n/(2*n+1);
    end
    LHS = abs(pi - series);
end
sp4c = k;

%% Problem 5: (a-c) 5 pts each, (d) 10 pts, 25 pts
omega = -10:-10:-70;
for n = 1:length(omega)
    [T{n},X{n},Y{n},Z{n},U{n},V{n},W{n}] = soccer(omega(n));
end
sp5a = T{1}(end);
sp5b = [X{1}(end) Y{1}(end) Z{1}(end)];
sp5c = [U{1}(end) V{1}(end) W{1}(end)];

figure(3); hold on;
cs = 'krbgmcy';
for n = 1:length(omega)
    plot3(X{n},Y{n},Z{n},cs(n),'LineWidth',2);    
    legend_string{n} = sprintf('omega = %2d',omega(n));
end
for n = 1:length(omega)
    plot3(X{n}(end),Y{n}(end),Z{n}(end),[cs(n) 'o'],'MarkerFaceColor',cs(n)); 
end
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
legend(legend_string,'Location','best');
title('Effect of spinning rate omega');
view(3); box on; grid on;
set(gca,'FontSize',14)
sp5d = 'See figure 3';
