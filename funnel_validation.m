% Get variables

model_settings = get_protocol_settings(37);


% Parameters funnel
r0 = sqrt(125/pi);      % Initial radius (mm)
k = (sqrt(440/pi)-r0)/10;       % Tapering constant


% Grid for z <= 0 (Linear)
% Spacing is 1/model_settings.stp
x_neg = model_settings.xmin : (1/model_settings.stp) : 0;

% Grid for z > 0 (Volume-Weighted)
% First, calculate total number of elements needed for the positive region
M_total = (model_settings.stp / (3*k*r0^2)) * ((r0 + k*model_settings.xmax)^3 - r0^3);
M_total = ceil(M_total); % Round up to nearest integer

% Create an array of element indices
m_array = 1:M_total;

% Calculate positions x_pos using inverse formula
x_pos = (1/k)*(((3*k*r0^2/model_settings.stp).*m_array+r0^3).^(1/3) - r0);

% Filter out any points that exceed xmax due to rounding
x_pos = x_pos(x_pos <= model_settings.xmax);

% Combine
x_grid = [x_neg, x_pos];



%% Generate signal
v_ecg_list = [5 10 20 30 40];
v_resp_list = [0.5 1 2 4 7];
v_sv_list = [0.1 0.2 0.3 0.4 0.5];


num_vels = 5;
N_runs = 1;
Signal_store_funnel = zeros(nfreqs, 125, 170);
Qvars_fit_funnel = zeros(nfreqs, num_vels^3, 6);

% Loop thorugh all frequencies
for frq = 1:nfreqs
    
    % Get frequecies and phase shifts
    fecg = f_ecg(frq);
    fresp = f_resp(frq);
    fsv = f_lf(frq);

    decg = d_ecg(frq);
    dresp = d_resp(frq);
    dsv = d_lf(frq);

    frqs_sim = [fecg fresp fsv];

    % Loop through velocities
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
         
            for j = 1:(model_settings.Nt + model_settings.ds) % Loop through  excitations
            
                if (and(i==1,j==1))
                    continue
                end
        
                tcurrent = model_settings.t_ex(j, i); % Current excitation
                [~, Icurrent_ecg] = min(abs(t_ecg - tcurrent));
    
        
                % Do ECG
                CI(j, i, 1) = (sin(theta_ecg(Icurrent_ecg)))/w_ecg; %w_ecg_(Icurrent_ecg); %
                SI(j, i, 1) = -(cos(theta_ecg(Icurrent_ecg)) - 1)/w_ecg; %w_ecg_(Icurrent_ecg); %
    
                CI(j, i, 2) = (sin(theta_resp(Icurrent_ecg)))/w_resp; %w_ecg_(Icurrent_ecg); %
                SI(j, i, 2) = -(cos(theta_resp(Icurrent_ecg)) - 1)/w_resp; %w_ecg_(Icurrent_ecg); %
    
                CI(j, i, 3) = (sin(theta_sv(Icurrent_ecg)))/w_sv; %w_ecg_(Icurrent_ecg); %
                SI(j, i, 3) = -(cos(theta_sv(Icurrent_ecg)) - 1)/w_sv; %w_ecg_(Icurrent_ecg); %
    
            end
        
        end
   

    
        tstart = model_settings.t_ex(1, 1); % First excitation
    

    
        f = @(xvars) get_signal_funnel(xvars, x_grid, model_settings.L, model_settings.Nt, model_settings.ds, model_settings.n_slices, model_settings.x0, model_settings.t_ex, ...
            frqs_sim, model_settings.T1, model_settings.T2, model_settings.TE, model_settings.theta, 1, r0, k);


        h = @(xvars) get_signal_pipe(xvars, model_settings.xmin, model_settings.xmax, model_settings.L, model_settings.Nt, model_settings.ds, model_settings.n_slices, model_settings.x0, model_settings.t_ex, CI, SI, ...
            model_settings.T1, model_settings.T2, model_settings.TE, model_settings.theta, 1, 0, model_settings.stp);

        
        % Variable vecg: Modulo arithmetic loops the values from 1 to 5
        vecg = v_ecg_list((mod(ids - 1, 5) + 1));
        
        % Variable vresp: Ceiling division creates blocks of 5
        vresp = v_resp_list(ceil((mod(ids - 1, 25) + 1) / 5));
        
        % Variable vsv: Ceiling division creates blocks of 25
        vsv = v_sv_list(ceil(ids / 25));
    
        A = [vecg decg vresp dresp vsv dsv];

        y = f(A);
        Signal_store_funnel(frq, ids, :) = y' ;





        obj_f = @(x) (sum((y - h(x)).^2)); 
        
        Ab = squeeze(Qvars_fit_funnel(frq, ids, :));
        obj_b = obj_f(Ab);


        i = 0;
        % Make fits
        while i < N_runs
                
            lb = [-50, -15*pi, -10, -15*pi, -1, -15*pi];
            ub = [50, 15*pi, 10, 15*pi, 1, 15*pi];
                
            try
                options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon, 'UseParallel', true);
            
        
                A = particleswarm(obj_f, 6, lb, ub, options);
                if isempty(A)
                    disp('Particleswarm returned empty solution');
                    A = [0 0 0 0 0 0];
                else
                    i = i + 1;
                end
    
                if obj_f(A) < obj_b
                    obj_b = obj_f(A);
                    Ab = A;
                    disp('Improved')
                end
    
            catch ME
                disp('Particleswarm failed');
                A = [0 0 0 0 0 0];
            end


        end

        
        Qvars_fit_funnel(frq, ids, :) = Ab;



        if false
            figure()
            plot(y), hold on, plot(h(Ab)), legend('True', 'Fit')
            title(['C: ' num2str(vecg) ' vs ' num2str(Ab(1)) ' R: ' num2str(vresp) ...
                ' vs ' num2str(Ab(3)) ' LF: ' num2str(vsv) ' vs ' num2str(Ab(5))])
            pause(1)
        end



    end

end



