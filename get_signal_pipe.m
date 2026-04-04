function S = get_signal_pipe(vars, xmin, xmax, L, Nt, ds, m, x0, t_ex, CI, SI, ...
    T1, T2, TE, theta, M0, Dpin, stp)

    
    N = round(stp*(xmax-xmin)) + 1; % Nr of points
    
    x = linspace(xmin, xmax, N); % Grid points
    Mz = M0*ones(1,N); % Magnetization
    S = zeros(Nt, 1); % Signal

    
    [~, I] = sort(t_ex(:));
    sz = [Nt + ds, m];
    [row_t,col_t] = ind2sub(sz,I);
    

    time_last = -inf;

    for i = 1:length(I)

        % Extract current time and slice
        timePoint = row_t(i);
        slice = col_t(i);
        time_current = t_ex(timePoint, slice);

        % Calculate displacement
        D = vars(1)*(CI(timePoint, slice, 1)*cos(vars(2)) - ...
            SI(timePoint, slice, 1)*sin(vars(2))) + ...
            vars(3)*(CI(timePoint, slice, 2)*cos(vars(4)) - ...
            SI(timePoint, slice, 2)*sin(vars(4))) + ...
            vars(5)*(CI(timePoint, slice, 3)*cos(vars(6)) - ...
            SI(timePoint, slice, 3)*sin(vars(6)));


        % Perform dispersion
        if Dpin > 0
            nconv = 10*stp+1;
            dt = sqrt(time_current - time_last);
            Dp = Dpin*100*stp^2/60; % points^2/s
            alpha = (nconv-1)/(2*sqrt(2*Dp)*dt);
            w = gausswin(nconv, alpha);
            w = w/sum(w);
            Mz = conv([Mz(floor(nconv/2):-1:1), Mz, Mz((length(Mz)):-1:(length(Mz)-floor(nconv/2)+1))], w,'valid');
        end

        
        % Relaxation
        deltaT = time_current - time_last;
        Mz = M0*(1 - exp(-(deltaT)/T1)) + Mz.*exp(-(deltaT)/T1);
        


        % Get positions to excite
        x_ex = x0(slice) - D;
        I_pos = and(x >= x_ex, x < (x_ex + L));


        % Extract signal
        if (slice == 1 && timePoint>ds)
            S(timePoint-ds) = mean(Mz(I_pos)*sind(theta))*exp(-TE/T2);
        end

        % Update parallel magnetization and time
        Mz(I_pos) = Mz(I_pos)*cosd(theta);
        time_last = time_current;
    end




    % Normalize signal
    S = S/mean(S);


end
