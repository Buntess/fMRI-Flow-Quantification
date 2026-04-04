% Get settings for model
model_settings = get_protocol_settings(37);


% Define frequency and phase
fecg = 0.97;
t_ecg = 0:0.01:600;
theta_ecg = 2*pi*fecg*t_ecg;


% Displacements will be calculated relative to the first excitation
tstart = model_settings.t_ex(1, 1);
[~, Istart_ecg] = min(abs(t_ecg - tstart));
w_ecg = 2*pi*fecg;
theta_start_ecg = theta_ecg(Istart_ecg);
theta_ecg = mod(theta_ecg - theta_start_ecg, 2*pi);



CI = zeros(model_settings.Nt + model_settings.ds, model_settings.n_slices, 1);
SI = zeros(model_settings.Nt + model_settings.ds, model_settings.n_slices, 1);

% Perform pre-integration
for i = 1:model_settings.n_slices % Loop through all slices
     
    for j = 1:(model_settings.Nt + model_settings.ds) % Loop through  excitations
    
        if (and(i==1,j==1))
            continue
        end

        tcurrent = model_settings.t_ex(j, i); % Current excitation
        [~, Icurrent_ecg] = min(abs(t_ecg - tcurrent));


        CI(j, i, 1) = (sin(theta_ecg(Icurrent_ecg)))/w_ecg;
        SI(j, i, 1) = -(cos(theta_ecg(Icurrent_ecg)) - 1)/w_ecg;

    end

end



f = @(x, fp) get_signal_profile(x, model_settings.xmin, model_settings.xmax, model_settings.L, ...
    model_settings.Nt, model_settings.ds, model_settings.n_slices, model_settings.x0, model_settings.t_ex,...
    CI, SI, model_settings.T1, model_settings.T2, model_settings.TE, model_settings.theta, 1, model_settings.stp, fp);


%%


v = [1, 4, 7, 10, 14, 20, 25];


v_est = zeros(3, length(v));
N = round(model_settings.stp*(model_settings.xmax-model_settings.xmin)) + 1; % Nr of points

for jj = 1:3

    % Define flow profile
    if jj == 1
        A = 1;
        B = 3;  
        r = sqrt(A^2 + (B^2 - A^2)*rand(1, N));
        pf = -(r.^2 - B^2 - (B^2 - A^2)/log(B/A).*log(r/B));
        pf = pf/mean(pf);
        pf = pf(randperm(length(pf)));
    elseif jj == 2
        pf = 2*linspace(0,1,N);
        pf = pf(randperm(length(pf)));
    else
        pf = ones(1,N);
    end


    % Loop thorugh all velocities
    for ids = 1:length(v) 

        % Signal to match
        vv = v(ids);
        S = f([vv, pi/2], pf);
        y = S'/mean(S);
        
        % Fit plug flow
        fp_hc = ones(1, N); 
        obj_f = @(x) (sum((y - f(x, fp_hc)').^2)); 
        
        obj_b = inf;
        i = 0;
        while i < 1

            % Make fit
            try
                lb = [-2*vv, -5*pi];
                ub = [2*vv, 5*pi];
                options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon, 'UseParallel', true);
            
                A = particleswarm(obj_f, 2, lb, ub, options);
                if isempty(A)
                    disp('Particleswarm returned empty solution');
                else
                    i = i + 1;
                end

            catch ME
                disp('Particleswarm failed');
                A = [0 0 0 0 0 0];
            end
        
            if obj_f(A) < obj_b
                obj_b = obj_f(A);
                Ab = A;
            end
            pause(1)
        end

        % Save velocity amplitude
        v_est(jj, ids) = abs(Ab(1));

    end
end


%% Plot results 

v = [1, 4, 7, 10, 14, 20, 25];
clmaps = [0, 64, 55; 33, 135, 224; 255 193 7; 252, 191, 0]/255;

figure(1)
clf()
hold on
plot(v, v_est(1, :), '-^', 'MarkerSize',7, 'LineWidth',2.5, 'Color', clmaps(1,:), 'MarkerFaceColor', clmaps(1,:))
plot(v, v_est(2, :), 'o-', 'MarkerSize',7, 'LineWidth',2.5, 'Color', clmaps(2,:), 'MarkerFaceColor', clmaps(2,:))
plot(v, v_est(3, :), '-square', 'MarkerSize',7, 'LineWidth',2.5, 'Color', clmaps(3,:), 'MarkerFaceColor', clmaps(3,:))
ax = gca();
ax.FontSize = fs;
xlabel('True mean velocity (mm/s)', 'FontSize', font_size)
ylabel('Fitted mean velocity (mm/s)', 'FontSize', font_size)
legend('Annulus', 'Cylinder', 'Plug flow', 'Location', 'northwest', 'FontSize', fs)

axis([0 27.2 0 27.2])

rA = corr(v', v_est(1, :)')
rC = corr(v', v_est(2, :)')



