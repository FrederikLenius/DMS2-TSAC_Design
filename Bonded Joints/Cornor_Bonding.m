clear; clc; close all 

%% Parameters

% Adheasive strength [MPa]
% sigma_allow = 36; % Tensile strength from other glue 

tau_allow = 12*0.8; 
sigma_allow = tau_allow; % Shear strength used 

%% Bonded area database [m^2]

area(1).PanelNames = 'Bottom left';
area(1).A_norm_x = 54.06;    % [mm^2]
area(1).A_norm_y = 841.28;   % [mm^2]
area(1).A_norm_z = 995.06;   % [mm^2]
% Normal directions to ignore because they are finger side faces
area(1).IgnoreNormal = ["x"];

area(2).PanelNames = 'Firewall bottom';
area(2).A_norm_x = 0;      % [mm^2]
area(2).A_norm_y = 1622.95;      % [mm^2]
area(2).A_norm_z = 0;      % [mm^2]
% Normal directions to ignore because they are finger side faces
area(2).IgnoreNormal = [];


area(3).PanelNames = 'Firewall left';
area(3).A_norm_x = 0;      % [mm^2]
area(3).A_norm_y = 0;      % [mm^2]
area(3).A_norm_z = 667.8;      % [mm^2]
% Normal directions to ignore because they are finger side faces
area(3).IgnoreNormal = [];

area(4).PanelNames = 'Rear left';
area(4).A_norm_x = 667.80;      % [mm^2]
area(4).A_norm_y = 50.56;      % [mm^2]
area(4).A_norm_z = 763.20;      % [mm^2]
% Normal directions to ignore because they are finger side faces
area(4).IgnoreNormal = ["y"];


% Load parameters
data = readtable("AnsysGlueForces.csv");
timeMap = ["Side"; "Backward"; "Forward"; "Up"; "Down"];
data.Time = timeMap(data.Time);
data.Properties.VariableNames("Time") = "LoadCase";

% Convert force columns to numeric [N]
data.X = str2double(strrep(string(data.X), ',', '.'));
data.Y = str2double(strrep(string(data.Y), ',', '.'));
data.Z = str2double(strrep(string(data.Z), ',', '.'));

% Convert timestep too, if needed
data.Timestep = str2double(strrep(string(data.Step), ',', '.'));
%% Preallocate result table
Results = table();

%% Loop through all defined bonded edges
for i = 1:numel(area)

    panelName = string(area(i).PanelNames);

    % Find all force rows belonging to this edge/Name
    idx = strcmpi(string(data.Name), panelName);

    if ~any(idx)
        warning('No force data found for panel/Name: %s', panelName);
        continue
    end

    % Extract force rows for this panel
    panelData = data(idx, :);

    % Force components [N]
    Fx = panelData.X;
    Fy = panelData.Y;
    Fz = panelData.Z;

    A_x = area(i).A_norm_x;
    A_y = area(i).A_norm_y;
    A_z = area(i).A_norm_z;   

    %% Normal-check areas [mm^2]
    A_normal_x = A_x;
    A_normal_y = A_y;
    A_normal_z = A_z;
    
    %% Ignore selected normal directions
    if isfield(area, 'IgnoreNormal') && ~isempty(area(i).IgnoreNormal)
    
        ignoreDirs = lower(string(area(i).IgnoreNormal));
    
        if any(ignoreDirs == "x")
            A_normal_x = 0;
        end
    
        if any(ignoreDirs == "y")
            A_normal_y = 0;
        end
    
        if any(ignoreDirs == "z")
            A_normal_z = 0;
        end
    end
    

    % Effective shear areas [mm^2]
    % For X-load, Y- and Z-normal faces are parallel to the load
    A_shear_x = A_y + A_z;
    A_shear_y = A_x + A_z;
    A_shear_z = A_x + A_y;

    %% Initialize stresses and safety factors as invalid/unused
    nRows = height(panelData);
    
    sigma_x = NaN(nRows, 1);
    sigma_y = NaN(nRows, 1);
    sigma_z = NaN(nRows, 1);
    
    tau_x = NaN(nRows, 1);
    tau_y = NaN(nRows, 1);
    tau_z = NaN(nRows, 1);
    
    SF_normal_x = NaN(nRows, 1);
    SF_normal_y = NaN(nRows, 1);
    SF_normal_z = NaN(nRows, 1);
    
    SF_shear_x = NaN(nRows, 1);
    SF_shear_y = NaN(nRows, 1);
    SF_shear_z = NaN(nRows, 1);
    
    %% Normal stresses [MPa]
    % Only calculate if the relevant normal area exists
    if A_normal_x > 0
        sigma_x = Fx ./ A_normal_x;
        SF_normal_x = sigma_allow ./ abs(sigma_x);
    end
    
    if A_normal_y > 0
        sigma_y = Fy ./ A_normal_y;
        SF_normal_y = sigma_allow ./ abs(sigma_y);
    end
    
    if A_normal_z > 0
        sigma_z = Fz ./ A_normal_z;
        SF_normal_z = sigma_allow ./ abs(sigma_z);
    end
    
    %% Effective shear stresses [MPa]
    % Only calculate if the relevant shear area exists
    if A_shear_x > 0
        tau_x = Fx ./ A_shear_x;
        SF_shear_x = tau_allow ./ abs(tau_x);
    end
    
    if A_shear_y > 0
        tau_y = Fy ./ A_shear_y;
        SF_shear_y = tau_allow ./ abs(tau_y);
    end
    
    if A_shear_z > 0
        tau_z = Fz ./ A_shear_z;
        SF_shear_z = tau_allow ./ abs(tau_z);
    end
    
    %% Controlling safety factor per load direction
    % Ignore missing/invalid modes
    SF_x = min([SF_normal_x, SF_shear_x], [], 2, 'omitnan');
    SF_y = min([SF_normal_y, SF_shear_y], [], 2, 'omitnan');
    SF_z = min([SF_normal_z, SF_shear_z], [], 2, 'omitnan');
    
    %% Identify controlling direction/mode
    SF_all = [ ...
        SF_normal_x, SF_shear_x, ...
        SF_normal_y, SF_shear_y, ...
        SF_normal_z, SF_shear_z];
    
    ModeNames = [ ...
        "Normal X", "Shear X", ...
        "Normal Y", "Shear Y", ...
        "Normal Z", "Shear Z"];
    
    [SF_control, controlIdx] = min(SF_all, [], 2, 'omitnan');
    ControlMode = ModeNames(controlIdx).';
    
    %% Handle rows where no valid stress mode existed
    invalidRows = all(isnan(SF_all), 2);
    
    SF_control(invalidRows) = NaN;
    ControlMode(invalidRows) = "No valid area";

    %% Build result table for this panel
    panelResults = table();

    panelResults.Name    = panelData.Name;
    panelResults.LoadCase = panelData.LoadCase;
    panelResults.Timestep = panelData.Timestep;

    panelResults.Fx = Fx;
    panelResults.Fy = Fy;
    panelResults.Fz = Fz;

    panelResults.A_norm_x_mm2 = repmat(A_x, height(panelData), 1);
    panelResults.A_norm_y_mm2 = repmat(A_y, height(panelData), 1);
    panelResults.A_norm_z_mm2 = repmat(A_z, height(panelData), 1);

    panelResults.A_shear_x_mm2 = repmat(A_shear_x, height(panelData), 1);
    panelResults.A_shear_y_mm2 = repmat(A_shear_y, height(panelData), 1);
    panelResults.A_shear_z_mm2 = repmat(A_shear_z, height(panelData), 1);

    panelResults.sigma_x_MPa = sigma_x;
    panelResults.sigma_y_MPa = sigma_y;
    panelResults.sigma_z_MPa = sigma_z;

    panelResults.tau_x_MPa = tau_x;
    panelResults.tau_y_MPa = tau_y;
    panelResults.tau_z_MPa = tau_z;

    panelResults.SF_normal_x = SF_normal_x;
    panelResults.SF_normal_y = SF_normal_y;
    panelResults.SF_normal_z = SF_normal_z;

    panelResults.SF_shear_x = SF_shear_x;
    panelResults.SF_shear_y = SF_shear_y;
    panelResults.SF_shear_z = SF_shear_z;

    panelResults.SF_x = SF_x;
    panelResults.SF_y = SF_y;
    panelResults.SF_z = SF_z;

    panelResults.SF_min = SF_control;
    panelResults.ControlMode = ControlMode;

    %% Append to full results table
    Results = [Results; panelResults];
end

%% Sort by lowest safety factor
Results = sortrows(Results, "SF_min", "ascend");

%% Export full results
outputFile = 'adhesive_safety_factors.csv';
writetable(Results, outputFile);

fprintf('\nResults exported to: %s\n', outputFile);

%% Print total bonded glue area

totalGlueArea_mm2 = 0;

fprintf('\nBonded glue area summary:\n');

for i = 1:numel(area)

    A_vals_mm2 = [area(i).A_norm_x, area(i).A_norm_y, area(i).A_norm_z];
    A_vals_mm2(isnan(A_vals_mm2)) = 0;

    A_edge_mm2 = sum(A_vals_mm2);
    totalGlueArea_mm2 = totalGlueArea_mm2 + A_edge_mm2;

    fprintf('  %-25s : %.2f mm^2\n', ...
        area(i).PanelNames, ...
        A_edge_mm2);
end

fprintf('---------------------------------------------\n');
fprintf('  %-25s : %.2f mm^2\n\n', ...
    'Total bonded area', ...
    totalGlueArea_mm2);