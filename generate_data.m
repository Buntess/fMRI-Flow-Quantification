%% Get variables

std_noise = 0.1; % Random noise to signal
mu_cardiac = 72/60; % Hz
std_cardiac = 0.1; % Hz
mu_resp = 0.28; % Hz
std_resp = 0.05; % Hz
mu_lf = 0.071; % Hz
std_lf = 0.02; % Hz


% Make cardiac, respiratory and low frequency frequencies and phase shifts
rng(100)
nfreqs = 3;

f_ecg = mu_cardiac + std_cardiac*randn(nfreqs, 1)
f_resp = mu_resp + std_resp*randn(nfreqs, 1)
f_lf = mu_lf + std_lf*randn(nfreqs, 1)


d_ecg = 2*pi*rand(nfreqs, 1)
d_resp = 2*pi*rand(nfreqs, 1)
d_lf = 2*pi*rand(nfreqs, 1)



% Font size for plots
fs = 18;
font_size = 26;

v_ecg_list = [5 10 20 30 40];
v_resp_list = [0.5 1 2 4 7];
v_sv_list = [0.1 0.2 0.3 0.4 0.5];


%% Create Signal_store (without noise)

model_settings = get_protocol_settings();
Dp = 6.0; % Dispersion coefficient

num_vels = 5;
v_ecg_list = [5 10 20 30 40];
v_resp_list = [0.5 1 2 4 7];
v_sv_list = [0.1 0.2 0.3 0.4 0.5];


Signal_store = zeros(nfreqs, num_vels^3, model_settings.Nt);

for frq = 1:nfreqs
 
    fecg = f_ecg(frq);
    fresp = f_resp(frq);
    fsv = f_lf(frq);

    decg = d_ecg(frq);
    dresp = d_resp(frq);
    dsv = d_lf(frq);
    ids = 0;
    for cvel = 1:(num_vels^3)
        ids = ids + 1
   

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
    

                CI(j, i, 1) = (sin(theta_ecg(Icurrent_ecg)))/w_ecg; %w_ecg_(Icurrent_ecg); %
                SI(j, i, 1) = -(cos(theta_ecg(Icurrent_ecg)) - 1)/w_ecg; %w_ecg_(Icurrent_ecg); %
    
                CI(j, i, 2) = (sin(theta_resp(Icurrent_ecg)))/w_resp; %w_ecg_(Icurrent_ecg); %
                SI(j, i, 2) = -(cos(theta_resp(Icurrent_ecg)) - 1)/w_resp; %w_ecg_(Icurrent_ecg); %
    
                CI(j, i, 3) = (sin(theta_sv(Icurrent_ecg)))/w_sv; %w_ecg_(Icurrent_ecg); %
                SI(j, i, 3) = -(cos(theta_sv(Icurrent_ecg)) - 1)/w_sv; %w_ecg_(Icurrent_ecg); %
    
            end
        
        end
        
        

    
        f = @(x) get_signal_pipe(x, model_settings.xmin, model_settings.xmax, model_settings.L,...
            model_settings.Nt, model_settings.ds, model_settings.n_slices, model_settings.x0, ...
            model_settings.t_ex, CI, SI, model_settings.T1, model_settings.T2, model_settings.TE, ...
            model_settings.theta, 1, Dp, model_settings.stp);
        
        % Variable vecg: Modulo arithmetic loops the values from 1 to 5
        vecg = v_ecg_list((mod(ids - 1, 5) + 1));
        
        % Variable vresp: Ceiling division creates blocks of 5
        vresp = v_resp_list(ceil((mod(ids - 1, 25) + 1) / 5));
        
        % Variable vsv: Ceiling division creates blocks of 25
        vsv = v_sv_list(ceil(ids / 25));
    
        A = [vecg decg vresp dresp vsv dsv];

        y = f(A);
        Signal_store(frq, ids, :) = y';
    


    end

end