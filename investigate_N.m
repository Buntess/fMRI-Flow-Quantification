% Get variables 

model_settings = get_protocol_settings(37);
Dp = 6.0; % Dispersion coefficient





% Calculate signal – with different spatial resolution

num_vels = 5;
v_ecg_list = [5 10 20 30 40];
v_resp_list = [0.5 1 2 4 7];
v_sv_list = [0.1 0.2 0.3 0.4 0.5];


Nsteps_list = [100 50 20 10 5 1];



%% Calculate signals

nNsteps = length(Nsteps_list);
Signal_spatial_res = zeros(nNsteps, num_vels^3, model_settings.Nt);


for nstp = 1:nNsteps


    nnn = Nsteps_list(nstp);

    % Get frequencies and phase shifts
    fecg = f_ecg(1);
    fresp = f_resp(1);
    fsv = f_lf(1);

    decg = d_ecg(1);
    dresp = d_resp(1);
    dsv = d_lf(1);


    % Loop thorugh velocities
    for ids = 1:(num_vels^3)
   

        CI = zeros(model_settings.Nt + model_settings.ds, model_settings.n_slices, 3);
        SI = zeros(model_settings.Nt + model_settings.ds, model_settings.n_slices, 3);
        
        
        % Introduce cardiac, respiratory and low frequency phases
        t_ecg = 0:0.01:600;
        theta_ecg = 2*pi*fecg*t_ecg;

        t_resp = 0:0.01:600;
        theta_resp = 2*pi*fresp*t_resp;

        t_sv = 0:0.01:600;
        theta_sv = 2*pi*fsv*t_sv;
    
        
    
        % Displacements will be calculated relative to the first excitation
        tstart = model_settings.t_ex(1, 1); % First excitation
    
        % Cardiac
        [~, Istart_ecg] = min(abs(t_ecg - tstart));
        w_ecg = 2*pi*fecg;
        
        theta_start_ecg = theta_ecg(Istart_ecg);
        theta_ecg = mod(theta_ecg - theta_start_ecg, 2*pi);
    
        % Respiratory
        [~, Istart_resp] = min(abs(t_resp - tstart));
        w_resp = 2*pi*fresp;
        
        theta_start_resp = theta_resp(Istart_resp);
        theta_resp = mod(theta_resp - theta_start_resp, 2*pi);
    

        % Low frequency
        [~, Istart_sv] = min(abs(t_sv - tstart));
        w_sv = 2*pi*fsv;
        
        theta_start_sv = theta_sv(Istart_sv);
        theta_sv = mod(theta_sv - theta_start_sv, 2*pi);
 
    
        % Perform pre-integration of each component
        for i = 1:model_settings.n_slices % Loop through all slices
         
            for j = 1:(model_settings.Nt + model_settings.ds) % Loop through excitations
            
                if (and(i==1,j==1))
                    continue
                end
        
                tcurrent = model_settings.t_ex(j, i); % Current excitation
                [~, Icurrent_ecg] = min(abs(t_ecg - tcurrent));
    
        
                % Do ECG
                CI(j, i, 1) = (sin(theta_ecg(Icurrent_ecg)))/w_ecg; 
                SI(j, i, 1) = -(cos(theta_ecg(Icurrent_ecg)) - 1)/w_ecg;
    
                CI(j, i, 2) = (sin(theta_resp(Icurrent_ecg)))/w_resp;
                SI(j, i, 2) = -(cos(theta_resp(Icurrent_ecg)) - 1)/w_resp;
    
                CI(j, i, 3) = (sin(theta_sv(Icurrent_ecg)))/w_sv;
                SI(j, i, 3) = -(cos(theta_sv(Icurrent_ecg)) - 1)/w_sv;
    
            end
        
        end
        
        

    
        f = @(x) get_signal_pipe(x, model_settings.xmin, model_settings.xmax, model_settings.L, model_settings.Nt,...
            model_settings.ds, model_settings.n_slices, model_settings.x0, model_settings.t_ex, CI, SI, ...
            model_settings.T1, model_settings.T2, model_settings.TE, model_settings.theta, 1, Dp, nnn);
        

        % Variable vecg: Modulo arithmetic loops the values from 1 to 5
        vecg = v_ecg_list((mod(ids - 1, 5) + 1));
        
        % Variable vresp: Ceiling division creates blocks of 5
        vresp = v_resp_list(ceil((mod(ids - 1, 25) + 1) / 5));
        
        % Variable vsv: Ceiling division creates blocks of 25
        vsv = v_sv_list(ceil(ids / 25));
    
        A = [vecg decg vresp dresp vsv dsv];

        y = f(A); 
        Signal_spatial_res(nstp, ids, :) = y';




    end

end

 










%% Results with plots
% =========================================================================
% Relative Error and 95% Confidence Interval Calculation (Across n_sims)
% =========================================================================

Nsteps_list = [100 50 20 10 5 1];
% --- Dynamically read dimensions from your simulated data ---
num_I = size(Signal_spatial_res, 1);    
num_sims = size(Signal_spatial_res, 2);  
num_points = size(Signal_spatial_res, 3); 

% --- Calculate Relative Error ---
% Preallocate output array for speed
rel_error = zeros(num_I, num_sims);

disp('Calculating relative error...');
for I = 1:num_I
    for K = 1:num_sims
        % Extract the timecourse/points as column vectors for a specific simulation K
        x = squeeze(Signal_spatial_res(I, K, :)); 
        y = squeeze(Signal_spatial_res(1, K, :));
        
        % Calculate relative error: ||x - y||_2 / ||y||_2
        rel_error(I, K) = 100 * norm(x - y) / norm(y);
    end
end
disp('Calculation complete.');

% --- Process Data for Plotting ---
% Calculate mean and standard deviation of relative error across simulations
rel_err_avg = mean(rel_error, 2);


% --- Plotting ---
figure();
x_axis = 1:num_I;

% Create the error bar plot
plot(x_axis, rel_err_avg(num_I:-1:1), '-o', 'LineWidth', 2.5, 'MarkerSize', 8,...
    'MarkerFaceColor', '#D95319', 'Color', '#D95319');

ax = gca();
ax.FontSize = fs;
% Formatting the plot
xlabel('Resolution Steps (Elements per mm)', 'FontSize', font_size);
ylabel('Average Relative Error in %', 'FontSize', font_size);
xlim([0.5, num_I + 0.5]);


% Add a reference line at perfect agreement (Relative Error = 0)
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
ylim([-0.5 3.0])

% Dynamically create tick labels if Nsteps_list is available
if exist('Nsteps_list', 'var') && length(Nsteps_list) >= num_I
    labels = arrayfun(@num2str, Nsteps_list(num_I:-1:1), 'UniformOutput', false);
    xticks(1:num_I);
    xticklabels(labels);
else
    xticks(1:num_I);
end


