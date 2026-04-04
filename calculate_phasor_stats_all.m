
function [A_avg, B_avg, std_A, std_B, X, Y] = calculate_phasor_stats_all(A, B)
    % Calculates the mean and standard deviation of amplitudes (A) and 
    % phases (B) and returns the individual Cartesian coordinates.

    
    % Ensure inputs are column vectors
    A = A(:); 
    B = B(:);
    
    % 1. Map to Cartesian Coordinates (The "All Values")
    X = A .* cos(B);
    Y = A .* sin(B);
    
    % 2. Calculate Sample Statistics
    mean_X = mean(X);
    mean_Y = mean(Y);
    
    var_X = var(X);
    var_Y = var(Y);
    
    % Covariance matrix
    cov_matrix = cov(X, Y);
    cov_XY = cov_matrix(1, 2);
    
    % Calculate average amplitude and phase from the centroid
    A_avg = sqrt(mean_X^2 + mean_Y^2);
    B_avg = atan2(mean_Y, mean_X);
    
    % 3. Apply the Delta Method for Variances
    cos_B = cos(B_avg);
    sin_B = sin(B_avg);
    
    % Approximation for Amplitude Variance
    var_A = (cos_B^2 * var_X) + (sin_B^2 * var_Y) + (2 * cos_B * sin_B * cov_XY);
    
    % Approximation for Phase Variance:
    var_B = ((sin_B^2 / A_avg^2) * var_X) + ...
            ((cos_B^2 / A_avg^2) * var_Y) - ...
            (2 * (sin_B * cos_B / A_avg^2) * cov_XY);
            
    % 4. Extract Standard Deviations
    std_A = sqrt(var_A);
    std_B = sqrt(var_B);
end

