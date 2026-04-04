
function S = get_signal_profile(vars, xmin, xmax, L, Nt, ds, m, x0, t_ex, CI, SI, ...
    T1, T2, TE, theta, M0, stp, v_dis)


    N = round(stp*(xmax-xmin)) + 1; % Nr of points

    xstart = linspace(xmin, xmax, N); % Grid points
    y = v_dis; % Flow profile


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
            SI(timePoint, slice, 1)*sin(vars(2)));

        % Find excitation region and update positions
        x_ex = x0(slice);
        x = xstart + D*y;
        I_pos = and(x >= x_ex, x < (x_ex + L));




        % Relaxation from last excitation
        deltaT = time_current - time_last;
        Mz = M0*(1 - exp(-(deltaT)/T1)) + Mz.*exp(-(deltaT)/T1);
        
        % Get signal
        if (slice == 1 && timePoint>ds)
            if sum(I_pos)>0
                S(timePoint-ds) = mean(Mz(I_pos)*sind(theta))*exp(-TE/T2);
            else
                S(timePoint-ds) = 0;
            end
        end

        % Uptade parallel magnetization and time
        Mz(I_pos) = Mz(I_pos)*cosd(theta);
        time_last = time_current;
    end

    % Normalize signal
    S = S/mean(S);


end
