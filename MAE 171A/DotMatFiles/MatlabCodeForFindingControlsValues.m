files_Mass1Still = {'Mass1Still1.mat', 'Mass1Still2.mat', 'Mass1Still3.mat', 'Mass1Still4.mat', 'Mass1Still5.mat'};
files_Mass2Still = {'Mass2Still1.mat', 'Mass2Still2.mat', 'Mass2Still3.mat', 'Mass2Still4.mat', 'Mass2Still5.mat'};
files_TwoDOF = {'TwoDOF1.mat', 'TwoDOF2.mat', 'TwoDOF3.mat', 'TwoDOF4.mat', 'TwoDOF5.mat'};

K1 = [];
K1plusK2 = [];
K2 = [];
M1 = [];
M2 = [];
D1 = [];
D2 = [];

F = 0.5;
n = 1;

for i = 1:5
    % Load TwoDOF data
    load(files_TwoDOF{i});
    Steady_State_disp = enc1(7000);
    K1(i) = F / Steady_State_disp;
end

k1avg = mean(K1);

for i = 1:5
    % Load Mass2Still data
    load(files_Mass2Still{i});
    Steady_State_disp = enc1(7000);
    K1plusK2(i) = F / Steady_State_disp;
end

K2 = K1plusK2 - K1;

for i = 1:5
    % Load Mass2Still data for M1
    load(files_Mass2Still{i});
    pks = findpeaks(enc1);
    om1 = find(enc1 == pks(5), 1);
    om2 = find(enc1 == pks(6), 1);
    t1 = om2 / 1000;
    t0 = om1 / 1000;
    wd = (2 * pi * n) / (t1 - t0);
    betaWn = 1 / (t1 - t0) * log((pks(1) - enc1(7000)) / (pks(2) - enc1(7000)));
    Wn = sqrt(wd^2 + betaWn^2);
    M1(i) = K1plusK2(i) / (Wn^2);
end

for i = 1:5
    % Load Mass1Still data for M2
    load(files_Mass1Still{i});
    pks = findpeaks(enc2);
    pks = pks(pks > 2500);
    om1 = find(enc2 == pks(1), 1);
    om2 = find(enc2 == pks(2), 1);
    t1 = om2 / 1000;
    t0 = om1 / 1000;
    wd = (2 * pi * n) / (t1 - t0);
    betaWn = 1 / (t1 - t0) * log((pks(1) - enc2(7000)) / (pks(2) - enc2(7000)));
    Wn = sqrt(wd^2 + betaWn^2);
    M2(i) = K2(i) / (Wn^2);
end

for i = 1:5
    % Load Mass2Still data for D1
    load(files_Mass2Still{i});
    pks = findpeaks(enc1);
    pks = pks(pks > 1750);
    t0 = find(enc1 == pks(1), 1) / 1000;
    t1 = find(enc1 == pks(2), 1) / 1000;
    wd = 2 * pi / (t1 - t0);
    betaWn = 1 / (t1 - t0) * log((pks(1) - enc1(7000)) / (pks(2) - enc1(7000)));
    Wn = sqrt(wd^2 + betaWn^2);
    beta = betaWn / Wn;
    D1(i) = K1plusK2(i) * 2 * beta / Wn;
end

for i = 1:5
    % Load Mass1Still data for D2
    load(files_Mass1Still{i});
    pks = findpeaks(enc2);
    pks = pks(pks > 1750);
    t0 = find(enc2 == pks(1), 1) / 1000;
    t1 = find(enc2 == pks(2), 1) / 1000;
    wd = 2 * pi / (t1 - t0);
    betaWn = 1 / (t1 - t0) * log((pks(1) - enc2(7000)) / (pks(2) - enc2(7000)));
    Wn = sqrt(wd^2 + betaWn^2);
    beta = betaWn / Wn;
    D2(i) = K2(i) * 2 * beta / Wn;
end

% Calculate averages
AvgK1 = mean(K1);
AvgK2 = mean(K2);
AvgM1 = mean(M1);
AvgM2 = mean(M2);
AvgD1 = mean(D1);
AvgD2 = mean(D2);

% Plot enc1 from the last loaded dataset
plot(enc1);
