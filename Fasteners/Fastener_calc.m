clear; clc; close all

%% Parameters
n_bolt  = 6;     % [-] Number of bolts

d_inner = 4.77;  % [mm] Inner bolt diameter
d_nom   = 6;     % [mm] Nominal bolt diameter
sig_UTS = 800;   % [MPa] Bolt UTS 
grf     = 0.8;   % [-] Bolt grade factor
A       = 20.1;  % [mm^2] Bolt stress area
T_lim   = 10.5;  % [Nm] Pre tension torque limit
K       = 0.125; % [-] Torque factor for dry bolts

m_TSAC  = 65;    % [kg] Total TSAC mass

% Friction coefficient
%mu = 0.5;   % [-] Steel-Steel
%mu = 0.47;  % [-] Alu-Steel 
%mu = 0.21;  % [-] Steel-Graphite
%mu = 1.35;  % [-] Alu-Alu
%mu = 0.18;  % [-] Alu-Rope
%mu = 0.14;  % [-] Steel-Carbon
mu = 0.33;  % [-] Kevlar-Steel, Figure 15 in https://link.springer.com/content/pdf/10.1007/BF01191969.pdf %
%mu = 0.20;  % [-] Polyethylene-steel 



V_TSAC = m_TSAC * 9.82 * 40; % [N] Total shear load of TSAC at 40g
sig_y = sig_UTS * grf;       % [MPa] Bolt yield strength

%% M6 chassis bolt pure shear and pre tension safety factor

% Pure shear
tau = V_TSAC / (A * n_bolt); % [MPa] Shear stress in each bolt
sig_vm = tau * sqrt(3);      % [MPa] Equivalent von Mieses stress

sf_y_shear = sig_y / sig_vm;       % [-] Safety factor for M6 bolt to yield
sf_UTS_shear = sig_UTS / sig_vm;   % [-] Safety factor for M6 bolt to UTS

% Pre tension
F_fric = V_TSAC/n_bolt;           % [F] Friction force per bolt
F_pre = V_TSAC/(mu*n_bolt);       % [F] Normal force, or pre tension force
F_pre_lim = 0.7 * sig_y * A;      % [F] Pre tension force limit

sig_pre = F_pre/A;                % [MPa] Pre tension stress
T_pre = K * F_pre * (d_nom*1e-3); % [Nm] Pre tension torque
sf_y_pre = sig_y / sig_pre;       % [-] Safety factor for M6 bolt to yield
sf_UTS_pre = sig_UTS / sig_pre;   % [-] Safety factor for M6 bolt to UTS
sf_Flim = F_pre_lim/F_pre;        % [-] Safety factor for M6 bolt force limit
sf_Tlim = T_lim/T_pre;            % [-] Safety factor for M6 bolt torque limit

%% Output
fprintf('--- Bolt Load Analysis (M6) ---\n')
fprintf('\nBolt paramaters ---------------\n')
fprintf('Bolt grade:    %6.1f\n', (sig_UTS/100+grf))
fprintf('Yield strength: %6.2f MPa\n', sig_y)
fprintf('UTS:            %6.2f MPa\n', sig_UTS)
fprintf('\n')

fprintf('Applied loads -----------------\n')
fprintf('Total shear force: %6.2f N\n',V_TSAC)
fprintf('Bolt shear force:  %6.2f N\n',V_TSAC/n_bolt)
fprintf('\n')

fprintf('\nPure shear ------------------\n')
fprintf('Shear stress per bolt: %6.2f MPa\n', tau)
fprintf('Von Mises stress:      %6.2f MPa\n', sig_vm)

fprintf('Safety factor (Yield, pure shear): %6.2f\n', sf_y_shear)
fprintf('Safety factor (UTS, pure shear):   %6.2f\n', sf_UTS_shear)
fprintf('\n')

fprintf('\nPre tension ------------------\n')
fprintf('Required friction force:  %6.2f N\n', F_fric)
fprintf('Pre tension force:        %6.2f N\n', F_pre)
fprintf('Pre tension force limit:  %6.2f N\n', F_pre_lim)
fprintf('\n')
fprintf('Pre tension stress:       %6.2f MPa\n', sig_pre)
fprintf('Pre tension torque:       %6.2f Nm\n', T_pre)
fprintf('Pre tension torque limit: %6.2f Nm\n', T_lim)
fprintf('\n')

fprintf('Safety factor (Yield, pre tension): %6.2f\n', sf_y_pre)
fprintf('Safety factor (UTS, pre tension):   %6.2f\n', sf_UTS_pre)
fprintf('Safety factor (Pre tension Force):  %6.2f\n', sf_Flim)
fprintf('Safety factor (Pre tension Torque): %6.2f\n', sf_Tlim)
fprintf('\n')
