function S = get_signal_funnel(vars, x, L, Nt, ds, num_slices, x0, t_ex, freqs, ...
    T1, T2, TE, theta, M0, r0, k)

    N = length(x);
    S = zeros(Nt, 1);

    [~, I] = sort(t_ex(:));
    sz = [Nt + ds, num_slices];
    [row_t,col_t] = ind2sub(sz,I);


    % Funnel geometry
    v_amp_scale = @(x) (r0./(r0 + k*x)).^2 .* (x>0) + (x<=0);


    % Initialize parallel magentizations and initial time
    Mz = M0 * ones(1, N); 
    tstart = t_ex(row_t(1), col_t(1));
    time_last = tstart; 
    
    for i = 1:length(I)

        % Extract current time and slice
        timePoint = row_t(i);
        slice = col_t(i);
        time_current = t_ex(timePoint, slice);
        
        v_time = @(t_val) vars(1)*cos(2*pi*freqs(1)*(t_val - tstart) + vars(2)) + ...
                          vars(3)*cos(2*pi*freqs(2)*(t_val - tstart) + vars(4)) + ...
                          vars(5)*cos(2*pi*freqs(3)*(t_val - tstart) + vars(6));
        
        v_func = @(pos, t_val) v_time(t_val) .* v_amp_scale(pos);
        

        % Get time step
        t_step = time_current - time_last; 
        
        % Update position via Runge–Kutta method
        if t_step > 0 
            t = time_last;
            k1 = v_func(x, t);
            k2 = v_func(x + (t_step/2)*k1, t + t_step/2);
            k3 = v_func(x + (t_step/2)*k2, t + t_step/2);
            k4 = v_func(x + t_step*k3, t + t_step);
            
            % Update position accurately
            x = x + (t_step/6) * (k1 + 2*k2 + 2*k3 + k4);
            
            % Apply T1 relaxation over the time step
            Mz = M0*(1 - exp(-t_step/T1)) + Mz.*exp(-t_step/T1);
        end
        
        % Find positions to be excited
        x_ex = x0(slice);
        I_pos = and(x >= x_ex, x < (x_ex + L));
        
        % Get signal
        if (slice == 1 && timePoint > ds)
            S(timePoint-ds) = mean(Mz(I_pos)*sind(theta))*exp(-TE/T2);
        end
        
        % Uptade parallel magnetization and time
        Mz(I_pos) = Mz(I_pos)*cosd(theta);
        time_last = time_current;
    end


    % Normalize signal
    S = S/mean(S);


end
