function [AvgK1, AvgK2, AvgM1, AvgM2, AvgD1, AvgD2] = compute_parameters(files_Mass1Still, files_Mass2Still, files_TwoDOF, F, n)
    % Initialize variables
    K1 = [];
    K1plusK2 = [];
    K2 = [];
    M1 = [];
    M2 = [];
    D1 = [];
    D2 = [];
    
    % Calculate K1
    for i = 1:5
        % Load TwoDOF data
        load(files_TwoDOF{i});
        Steady_State_disp = enc1(7000);
        K1(i) = F / Steady_State_disp;
    end

    k1avg = mean(K1);

    % Calculate K2
    for i = 1:5
        % Load Mass2Still data
        load(files_Mass2Still{i});
        Steady_State_disp = enc1(7000);
        K1plusK2(i) = F / Steady_State_disp;
    end

    K2 = K1plusK2 - K1;

    % Calculate M1
    for i = 1:5
        % Load Mass2Still data for M1
        load(files_Mass2Still{i});
        pks = findpeaks(enc1);

        if i == 3
            om1 = find(enc1 == pks(5), 1);
            om2 = find(enc1 == pks(6), 1);
        else
            om1 = find(enc1 == pks(1), 1);
            om2 = find(enc1 == pks(n+1), 1);
        end
        t1 = om2 / 1000;
        t0 = om1 / 1000;
        wd = (2 * pi * n) / (t1 - t0);
        betaWn = (1 / (t1 - t0)) * log((pks(1) - enc1(7000)) / (pks(2) - enc1(7000)));
        Wn = sqrt(wd^2 + betaWn^2);
        M1(i) = K1plusK2(i) / (Wn^2);
    end

    % Calculate M2
    for i = 1:5
        % Load Mass1Still data for M2
        load(files_Mass1Still{i});
        pks = findpeaks(enc2);
        pks = pks(pks > 2500);
        om1 = find(enc2 == pks(1), 1);
        om2 = find(enc2 == pks(n+1), 1);
        t1 = om2 / 1000;
        t0 = om1 / 1000;
        wd = (2 * pi * n) / (t1 - t0);
        betaWn = 1 / (t1 - t0) * log((pks(1) - enc2(7000)) / (pks(2) - enc2(7000)));
        Wn = sqrt(wd^2 + betaWn^2);
        M2(i) = K2(i) / (Wn^2);
    end

    % Calculate D1
    for i = 1:5
        % Load Mass2Still data for D1
        load(files_Mass2Still{i});
        pks = findpeaks(enc1);
        pks = pks(pks > 1750);
        t0 = find(enc1 == pks(1), 1) / 1000;
        t1 = find(enc1 == pks(n+1), 1) / 1000;
        wd = 2 * pi/ (t1 - t0);
        betaWn = 1 / (t1 - t0) * log((pks(1) - enc1(7000)) / (pks(2) - enc1(7000)));
        Wn = sqrt(wd^2 + betaWn^2);
        beta = betaWn / Wn;
        D1(i) = K1plusK2(i) * 2 * beta / Wn;
    end

    % Calculate D2
    for i = 1:5
        % Load Mass1Still data for D2
        load(files_Mass1Still{i});
        pks = findpeaks(enc2);
        pks = pks(pks > 1750);
        t0 = find(enc2 == pks(1), 1) / 1000;
        t1 = find(enc2 == pks(n+1), 1) / 1000;
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
end
