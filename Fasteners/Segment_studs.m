clear; clc; close all

%% Force calculation

Ms = 7.5;            % [kg] Segment mass
g  = 9.82;           % [m/s^2] Gravitational acceleration
ns = 4;              % [-] Number of studs
SF = 1;              % [-] Safety factor

Av = 20 * g;         % [m/s^2] Vertical acceleration
Al = 40 * g;         % [m/s^2] Lateral acceleration

F = Av * Ms;         % [N] Total pullout force
V = Al * Ms;         % [N] Total shear force
Fs = F/ns;           % [N] Pullout force pr stud
Vs = V/ns;           % [N] Shear force pr stud


%% Epoxy glue database

epoxies(1).name   = 'MG 9200 FR*';
epoxies(1).shear  = 10e6;  % [Pa] Al - Al

epoxies(2).name   = 'Permabond TA4230';
epoxies(2).shear  = 12e6*0.8;  % [Pa] CF - Cf (80% strength at 60 degC)

epoxies(3).name   = '3M DP100FR*';
epoxies(3).shear  = 07e6;  % [Pa] Al - Al

epoxies(4).name   = 'Delo FR898*';
epoxies(4).shear  = 18e6;  % [Pa] Al - Al

epoxies(5).name   = 'SR 1126 / SD 8205';
epoxies(5).shear  = 34*0.5*1e6;  % [Pa] Shear strength from tensile strength, converted with Tresca / Mohr's circle.


%% Stud database
% Areas found in CAD on Bossard website

studs(1).name  = 'M1/B20-M5x30';
studs(1).A_top = (152+55+21)*1e-6;         % [m^2] Top bonded area
studs(1).A_bot = (214+21+7*6)*1e-6;        % [m^2] Bottom bonded area
studs(1).A_bea = (0.5*pi*1.2*(9.8 + 20.0))*1e-6; % [m^2] Inserted side profile area

studs(2).name  = 'M1/B23-M5x30';
studs(2).A_top = (213+47+18)*1e-6;         % [m^2]
studs(2).A_bot = (272+18+7*5)*1e-6;        % [m^2]
studs(2).A_bea = (0.5*pi*1.2*(9.8 + 23.0))*1e-6; % [m^2] Inserted side profile area

studs(3).name  = 'M1/B38A-M5x30';
studs(3).A_top = (725+55+35+6*11)*1e-6;    % [m^2]
studs(3).A_bot = (798+35+6*11)*1e-6;       % [m^2]
studs(3).A_bea = (0.5*pi*1.2*(9.8 + 38.0))*1e-6; % [m^2] Inserted side profile area

%% Strength calculation

npass = 0;
ncomb = 0;
results = [];
all_results = [];
stud_results = [];

for i = 1:length(studs)
    A_top  = studs(i).A_top;           % [m^2] Top bonded area
    A_bot  = studs(i).A_bot;           % [m^2] Bottom bonded area
    A_tot   = A_top + A_bot;            % [m^2] Total bonded area
    A_bea = studs(i).A_bea;            % [m^2] Inserted profile side area

    % Embedded insert
    sigma_req_emb  = Fs / A_tot;            % [Pa] Required pullout stress
    tau_req_emb    = Vs / A_tot;            % [Pa] Required shear stress
    % Surface bonded insert
    sigma_req_sur  = Fs / A_bot;            % [Pa] Required pullout stress
    tau_req_sur    = Vs / A_bot;            % [Pa] Required shear stress

    sigma_bea = Vs / A_bea;            % [Pa] Bearing stress 

    stud_results(i).name = studs(i).name;
    stud_results(i).sig_bea = sigma_bea;
    stud_results(i).Area_top = A_top;
    stud_results(i).Area_bot = A_bot;
    stud_results(i).Area_tot = A_tot;
    for j = 1:length(epoxies)
        sigma_allow = epoxies(j).shear / SF;   % [Pa] Allowable tensile stress
        tau_allow   = epoxies(j).shear / SF;     % [Pa] Allowable shear stress

        SF_pullout_emb = sigma_allow / sigma_req_emb;
        SF_shear_emb   = tau_allow / tau_req_emb;

        pullout_ok_emb = sigma_req_emb <= sigma_allow;
        shear_ok_emb   = tau_req_emb   <= tau_allow;

        SF_pullout_sur = sigma_allow / sigma_req_sur;
        SF_shear_sur   = tau_allow / tau_req_sur;

        pullout_ok_sur = sigma_req_sur <= sigma_allow;
        shear_ok_sur   = tau_req_sur   <= tau_allow;
        
        % Store all combinations
        ncomb = ncomb + 1;

        all_results(ncomb).Stud = studs(i).name;
        all_results(ncomb).Epoxy = epoxies(j).name;

        all_results(ncomb).SF_pullout_emb = SF_pullout_emb;
        all_results(ncomb).SF_shear_emb   = SF_shear_emb;
        all_results(ncomb).SF_min_emb     = min(SF_pullout_emb, SF_shear_emb);
        all_results(ncomb).Pass_emb       = pullout_ok_emb && shear_ok_emb;

        all_results(ncomb).SF_pullout_sur = SF_pullout_sur;
        all_results(ncomb).SF_shear_sur   = SF_shear_sur;
        all_results(ncomb).SF_min_sur     = min(SF_pullout_sur, SF_shear_sur);
        all_results(ncomb).Pass_sur       = pullout_ok_sur && shear_ok_sur;

        if pullout_ok_emb && shear_ok_emb
            npass = npass + 1;

            results(npass).stud = studs(i).name;
            results(npass).epoxy = epoxies(j).name;
            results(npass).SF_pullout_emb = SF_pullout_emb;
            results(npass).SF_shear_emb   = SF_shear_emb;
            results(npass).SF_pullout_sur = SF_pullout_sur;
            results(npass).SF_shear_sur   = SF_shear_sur;
        end
    end
end

for l = 1:length(epoxies)
    Area_min = SF * max(Fs / epoxies(l).shear, Vs / epoxies(l).shear);  % [m^2] Minimum required bonding area
    epoxy_results(l).name = epoxies(l).name;
    epoxy_results(l).Area_min = Area_min;
end

%% Console Output

fprintf('Force calculations:\n')
fprintf('Total pullout force: %4.2f N\n', F)
fprintf('Total shear force: %4.2f N\n', V)
fprintf('Pullout force pr. stud: %4.2f N\n', Fs)
fprintf('Shear force pr. stud: %4.2f N\n', Vs)
fprintf('\n')

fprintf('\n')
fprintf('Minimum required bonding area for each glue\n')
for j = 1:length(epoxy_results)
    fprintf('%s | Area = %3.2f mm^2\n', ...
        epoxy_results(j).name, ...
        epoxy_results(j).Area_min * 1e6)
end

fprintf('\n')
fprintf('Stud parameters\n')
for i = 1:length(stud_results)
    fprintf('%s | Total Area = %3.2f mm^2 | Bottom area = %3.2f mm^2 | Top area %3.2f mm^2 | Bearing stress %3.2f MPa\n', ...
        stud_results(i).name, ...
        stud_results(i).Area_tot * 1e6, ...
        stud_results(i).Area_bot * 1e6, ...
        stud_results(i).Area_top * 1e6, ...
        stud_results(i).sig_bea * 1e-6)
end

fprintf('\n')
fprintf('Embedded stud + epoxy combinations that pass load check:\n')

if npass == 0
    fprintf('No combinations pass the load check.\n')
else
    % Sort by lowest safety factor first
    % (most critical acceptable design first)

    min_SF = zeros(1, npass);

    for k = 1:npass
        min_SF(k) = min(results(k).SF_pullout_emb, results(k).SF_shear_emb);
    end

    [~, sort_idx] = sort(min_SF, 'ascend');
    results = results(sort_idx);

    for k = 1:npass
        fprintf('%s + %s | pullout SF = %3.2f | shear SF = %3.2f\n', ...
            results(k).stud, ...
            results(k).epoxy, ...
            results(k).SF_pullout_emb, ...
            results(k).SF_shear_emb ...
            )
    end
end

fprintf('\n')
fprintf('Surface bonded stud + epoxy combinations that pass load check:\n')

if npass == 0
    fprintf('No combinations pass the load check.\n')
else
    % Sort by lowest safety factor first
    % (most critical acceptable design first)

    min_SF = zeros(1, npass);

    for h = 1:npass
        min_SF(h) = min(results(h).SF_pullout_sur, results(h).SF_shear_sur);
    end

    [~, sort_idx] = sort(min_SF, 'ascend');
    results = results(sort_idx);

    for h = 1:npass
        fprintf('%s + %s | pullout SF = %3.2f | shear SF = %3.2f\n', ...
            results(h).stud, ...
            results(h).epoxy, ...
            results(h).SF_pullout_sur, ...
            results(h).SF_shear_sur ...
            )
    end
end

fprintf('\n')
fprintf('*Shear strength not from composite bonding\n')
fprintf('SR 1126 uses shear stress from tresca/mohr''s circle\n')

%% Export to Excel

filename = 'Stud bonding results.xlsx';

% --- Epoxy sheet ---

epoxy_names = {epoxies.name}';
shear_MPa   = [epoxies.shear]' * 1e-6;  % [MPa]

Area_min_mm2 = [epoxy_results.Area_min]' * 1e6; % [mm^2]

T_epoxy = table(epoxy_names, shear_MPa, Area_min_mm2, ...
    'VariableNames', {'Epoxy', 'Shear Strength [MPa]', 'Min Area [mm^2]'});

sheet = 'Epoxy';

writetable(T_epoxy, filename, 'Sheet', sheet);

% Write disclaimers below table
nrows = height(T_epoxy) + 3;

writecell({'*Shear strength not from composite bonding'}, filename, ...
    'Sheet', sheet, 'Range', ['A' num2str(nrows)]);

writecell({'SR 1126 uses shear stress from Tresca/Mohr''s circle'}, filename, ...
    'Sheet', sheet, 'Range', ['A' num2str(nrows+1)]);


% --- Studs sheet ---
stud_names = {stud_results.name}';

A_top_mm2 = [stud_results.Area_top]' * 1e6;
A_bot_mm2 = [stud_results.Area_bot]' * 1e6;
A_tot_mm2 = [stud_results.Area_tot]' * 1e6;

A_tot_all_mm2 = [stud_results.Area_tot]' * ns * 1e6;
A_bot_all_mm2 = [stud_results.Area_bot]' * ns * 1e6;

sig_bea_MPa = [stud_results.sig_bea]' * 1e-6;

T_studs = table(stud_names, A_tot_mm2, A_bot_mm2, A_top_mm2, ...
    A_tot_all_mm2, A_bot_all_mm2, sig_bea_MPa, ...
    'VariableNames', {'Stud', 'AreaTotal_mm2', 'AreaBottom_mm2', ...
    'AreaTop_mm2', 'AreaTotal_allStuds_mm2', ...
    'AreaBottom_allStuds_mm2', 'BearingStress_MPa'});

writetable(T_studs, filename, 'Sheet', 'Studs');

% --- All combinations sheet ---
T_all = struct2table(all_results);
T_all = sortrows(T_all, 'SF_min_emb', 'ascend');
T_all.Properties.VariableNames = { ...
    'Stud', ...
    'Epoxy', ...
    'SF Pullout (Embedded)', ...
    'SF Shear (Embedded)', ...
    'SF Min (Embedded)', ...
    'Pass (Embedded)', ...
    'SF Pullout (Surface)', ...
    'SF Shear (Surface)', ...
    'SF Min (Surface)', ...
    'Pass (Surface)'};
sheet = 'All combinations';
writetable(T_all, filename, 'Sheet', sheet);

% Write disclaimers below table
nrows = height(T_all) + 3;

writecell({'*Shear strength not from composite bonding'}, filename, ...
    'Sheet', sheet, 'Range', ['A' num2str(nrows)]);

writecell({'SR 1126 uses shear stress from Tresca/Mohr''s circle'}, filename, ...
    'Sheet', sheet, 'Range', ['A' num2str(nrows+1)]);

fprintf('Results exported to %s\n', filename);