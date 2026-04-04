% Get variables

model_settings = get_protocol_settings(10);
Dp = 6.0;


%% Calc CI and SI

N_runs = 1;
Signal_store_pipe = zeros(nfreqs, num_vels^3, model_settings.Nt);
Qvars_fit = zeros(nfreqs, num_vels^3, 6);

for frq = 1:nfreqs
    frq
    % Get frequencies and phase shifts
    fecg = f_ecg(frq);
    fresp = f_resp(frq);
    fsv = f_lf(frq);

    decg = d_ecg(frq);
    dresp = d_resp(frq);
    dsv = d_lf(frq);


    % Loop thorugh velocities
    for ids = 1:(num_vels^3)
        ids

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
        
       


        % Variable vecg: Modulo arithmetic loops the values from 1 to 5
        vecg = v_ecg_list((mod(ids - 1, 5) + 1));
        
        % Variable vresp: Ceiling division creates blocks of 5
        vresp = v_resp_list(ceil((mod(ids - 1, 25) + 1) / 5));
        
        % Variable vsv: Ceiling division creates blocks of 25
        vsv = v_sv_list(ceil(ids / 25));

        



        f = @(xvars) get_signal_pipe(xvars, model_settings.xmin, model_settings.xmax, model_settings.L, model_settings.Nt, model_settings.ds, model_settings.n_slices, model_settings.x0, model_settings.t_ex, CI, SI, ...
            model_settings.T1, model_settings.T2, model_settings.TE, model_settings.theta, 1, Dp, model_settings.stp);



        %Generate data with noise
        A = [vecg decg vresp dresp vsv dsv];
        y = f(A);
        y = y + std_noise*randn(model_settings.Nt, 1);
        y = y/mean(y);
        Signal_store_pipe(frq, ids, :) = y';



        obj_f = @(x) (sum((y - f(x)).^2)); 
        
        Ab = squeeze(Qvars_fit(frq, ids, :));
        obj_b = obj_f(Ab);



        i = 0;
        % Make fit
        while i < N_runs
                
            lb = [-50, -15*pi, -10, -15*pi, -1, -15*pi];
            ub = [50, 15*pi, 10, 15*pi, 1, 15*pi];
                
            try
                options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon, 'UseParallel', true);
            
        
                A = particleswarm(obj_f, 6, lb, ub, options);
                if isempty(A)
                    disp('Particleswarm returned empty solution');
                    % Här kan du lägga till logik för att hantera detta fall
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

        
        Qvars_fit(frq, ids, :) = Ab;

        if true
            figure()
            plot(y), hold on, plot(f(Ab)), legend('True', 'Fit')
            title(['C: ' num2str(vecg) ' vs ' num2str(Ab(1)) ' R: ' num2str(vresp) ...
                ' vs ' num2str(Ab(3)) ' LF: ' num2str(vsv) ' vs ' num2str(Ab(5))])
            pause(1)


        end


    end

end

