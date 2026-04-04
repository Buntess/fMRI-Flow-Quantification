%% Get indeces with matching velocities
Lindex = 1:125;

Ivecg1 = (mod(Lindex - 1, 5) + 1) == 1;
Ivecg2 = (mod(Lindex - 1, 5) + 1) == 2;
Ivecg3 = (mod(Lindex - 1, 5) + 1) == 3;
Ivecg4 = (mod(Lindex - 1, 5) + 1) == 4;
Ivecg5 = (mod(Lindex - 1, 5) + 1) == 5;

Ivresp1 = ceil((mod(Lindex - 1, 25) + 1) / 5) == 1;
Ivresp2 = ceil((mod(Lindex - 1, 25) + 1) / 5) == 2;
Ivresp3 = ceil((mod(Lindex - 1, 25) + 1) / 5) == 3;
Ivresp4 = ceil((mod(Lindex - 1, 25) + 1) / 5) == 4;
Ivresp5 = ceil((mod(Lindex - 1, 25) + 1) / 5) == 5;

Ivsv1 = ceil(Lindex / 25) == 1;
Ivsv2 = ceil(Lindex / 25) == 2;
Ivsv3 = ceil(Lindex / 25) == 3;
Ivsv4 = ceil(Lindex / 25) == 4;
Ivsv5 = ceil(Lindex / 25) == 5;




A1 = squeeze((Qvars_fit(1,:,:)));
A2 = squeeze((Qvars_fit(2,:,:)));
A3 = squeeze(Qvars_fit(3,:,:));





clr_map = [0.1059, 0.6196, 0.4667; 0.8510, 0.3725, 0.0078; 0.4588, 0.4392, 0.7020];


%%

[CA1, ~, std_A1, ~, PCA1, ~] = calculate_phasor_stats_all([A1(Ivecg1,1); A2(Ivecg1,1); A3(Ivecg1,1)],...
    [A1(Ivecg1,2)- d_ecg(1); A2(Ivecg1,2)- d_ecg(2); A3(Ivecg1,2)- d_ecg(3)]);
[CA2, ~, std_A2, ~, PCA2, ~] = calculate_phasor_stats_all([A1(Ivecg2,1); A2(Ivecg2,1); A3(Ivecg2,1)],...
    [A1(Ivecg2,2)- d_ecg(1); A2(Ivecg2,2)- d_ecg(2); A3(Ivecg2,2)- d_ecg(3)]);
[CA3, ~, std_A3, ~, PCA3, ~] = calculate_phasor_stats_all([A1(Ivecg3,1); A2(Ivecg3,1); A3(Ivecg3,1)],...
    [A1(Ivecg3,2)- d_ecg(1); A2(Ivecg3,2)- d_ecg(2); A3(Ivecg3,2)- d_ecg(3)]);
[CA4, ~, std_A4, ~, PCA4, ~] = calculate_phasor_stats_all([A1(Ivecg4,1); A2(Ivecg4,1); A3(Ivecg4,1)],...
    [A1(Ivecg4,2)- d_ecg(1); A2(Ivecg4,2)- d_ecg(2); A3(Ivecg4,2)- d_ecg(3)]);
[CA5, ~, std_A5, ~, PCA5, ~] = calculate_phasor_stats_all([A1(Ivecg5,1); A2(Ivecg5,1); A3(Ivecg5,1)],...
    [A1(Ivecg5,2)- d_ecg(1); A2(Ivecg5,2)- d_ecg(2); A3(Ivecg5,2)- d_ecg(3)]);

figure()
hold on 
axis([0 50 0 50])
plot([0 50], [0 50], '--k', 'LineWidth', 0.75)


% Get the means
y_means = [CA1, CA2, CA3, CA4, CA5];

% Calculate the 95% Confidence Interval for each group
y_ci = [1.96 * std_A1, 1.96 * std_A2, 1.96 * std_A3, 1.96 * std_A4, 1.96 * std_A5];

% Plot the line with vertical error bars
errorbar(v_ecg_list, y_means, y_ci, '.', 'MarkerSize',26, 'LineWidth',2.75, 'MarkerFaceColor', clr_map(2, :), 'Color', clr_map(2,:));

ax = gca();
ax.FontSize = fs;
xlabel('True velocity (mm/s)', 'FontSize',font_size)
ylabel('Simulated velocity (mm/s)', 'FontSize',font_size)




% Correlations
XC = [v_ecg_list(1)*ones(size(PCA1)); v_ecg_list(2)*ones(size(PCA1)); v_ecg_list(3)*ones(size(PCA1)); ...
    v_ecg_list(4)*ones(size(PCA1)); v_ecg_list(5)*ones(size(PCA1))];
YC = [PCA1; PCA2; PCA3; PCA4; PCA5];
[r_C, ~] = corr(XC, YC)





% Respiration
figure()
hold on 

axis([0 9 -1 9])
plot([0 9], [0 9], '--k', 'LineWidth', 0.75)

[CA1, ~, std_A1, ~, PRA1, ~] = calculate_phasor_stats_all([A1(Ivresp1,3); A2(Ivresp1,3); A3(Ivresp1,3)],...
    [A1(Ivresp1,4)- d_resp(1); A2(Ivresp1,4)- d_resp(2); A3(Ivresp1,4)- d_resp(3)]);
[CA2, ~, std_A2, ~, PRA2, ~] = calculate_phasor_stats_all([A1(Ivresp2,3); A2(Ivresp2,3); A3(Ivresp2,3)],...
    [A1(Ivresp2,4)- d_resp(1); A2(Ivresp2,4)- d_resp(2); A3(Ivresp2,4)- d_resp(3)]);
[CA3, ~, std_A3, ~, PRA3, ~] = calculate_phasor_stats_all([A1(Ivresp3,3); A2(Ivresp3,3); A3(Ivresp3,3)],...
    [A1(Ivresp3,4)- d_resp(1); A2(Ivresp3,4)- d_resp(2); A3(Ivresp3,4)- d_resp(3)]);
[CA4, ~, std_A4, ~, PRA4, ~] = calculate_phasor_stats_all([A1(Ivresp4,3); A2(Ivresp4,3); A3(Ivresp4,3)],...
    [A1(Ivresp4,4)- d_resp(1); A2(Ivresp4,4)- d_resp(2); A3(Ivresp4,4)- d_resp(3)]);
[CA5, ~, std_A5, ~, PRA5, ~] = calculate_phasor_stats_all([A1(Ivresp5,3); A2(Ivresp5,3); A3(Ivresp5,3)],...
    [A1(Ivresp5,4)- d_resp(1); A2(Ivresp5,4)- d_resp(2); A3(Ivresp5,4)- d_resp(3)]);


% Get the means 
y_means = [CA1, CA2, CA3, CA4, CA5];

% Calculate the 95% Confidence Interval for the group
y_ci = [1.96 * std_A1, 1.96 * std_A2, 1.96 * std_A3, 1.96 * std_A4, 1.96 * std_A5];

% Plot the line with vertical error bars
errorbar(v_resp_list, y_means, y_ci, '.', 'MarkerSize',26, 'LineWidth',2.75, 'MarkerFaceColor', clr_map(1, :), 'Color', clr_map(1,:));

ax = gca();
ax.FontSize = fs;
xlabel('True velocity (mm/s)', 'FontSize',font_size)
ylabel('Simulated velocity (mm/s)', 'FontSize',font_size)



% Correlations
XR = [v_resp_list(1)*ones(size(PRA1)); v_resp_list(2)*ones(size(PRA1)); v_resp_list(3)*ones(size(PRA1)); ...
    v_resp_list(4)*ones(size(PRA1)); v_resp_list(5)*ones(size(PRA1))];
YR = [PRA1; PRA2; PRA3; PRA4; PRA5];
[r_R, ~] = corr(XR, YR)








% Low frequency
figure()
hold on 
axis([0 0.6 -0.2 0.8])
plot([0 0.8], [0 0.8], '--k', 'LineWidth', 0.75)



[CA1, ~, std_A1, ~, PLFA1, ~] = calculate_phasor_stats_all([A1(Ivsv1,5); A2(Ivsv1,5); A3(Ivsv1,5)],...
    [A1(Ivsv1,6)- d_lf(1); A2(Ivsv1,6)- d_lf(2); A3(Ivsv1,6)- d_lf(3)]);
[CA2, ~, std_A2, ~, PLFA2, ~] = calculate_phasor_stats_all([A1(Ivsv2,5); A2(Ivsv2,5); A3(Ivsv2,5)],...
    [A1(Ivsv2,6)- d_lf(1); A2(Ivsv2,6)- d_lf(2); A3(Ivsv2,6)- d_lf(3)]);
[CA3, ~, std_A3, ~, PLFA3, ~] = calculate_phasor_stats_all([A1(Ivsv3,5); A2(Ivsv3,5); A3(Ivsv3,5)],...
    [A1(Ivsv3,6)- d_lf(1); A2(Ivsv3,6)- d_lf(2); A3(Ivsv3,6)- d_lf(3)]);
[CA4, ~, std_A4, ~, PLFA4, ~] = calculate_phasor_stats_all([A1(Ivsv4,5); A2(Ivsv4,5); A3(Ivsv4,5)],...
    [A1(Ivsv4,6)- d_lf(1); A2(Ivsv4,6)- d_lf(2); A3(Ivsv4,6)- d_lf(3)]);
[CA5, ~, std_A5, ~, PLFA5, ~] = calculate_phasor_stats_all([A1(Ivsv5,5); A2(Ivsv5,5); A3(Ivsv5,5)],...
    [A1(Ivsv5,6)- d_lf(1); A2(Ivsv5,6)- d_lf(2); A3(Ivsv5,6)- d_lf(3)]);

% Get the means
y_means = [CA1, CA2, CA3, CA4, CA5];


% Calculate the 95% Confidence Interval of each group
y_ci = [1.96 * std_A1, 1.96 * std_A2, 1.96 * std_A3, 1.96 * std_A4, 1.96 * std_A5];

% Plot the line with vertical error bars
errorbar(v_sv_list, y_means, y_ci, '.', 'MarkerSize',26, 'LineWidth',2.75, 'MarkerFaceColor', clr_map(3, :), 'Color', clr_map(3,:));

ax = gca();
ax.FontSize = fs;
xlabel('True velocity (mm/s)', 'FontSize',font_size)
ylabel('Simulated velocity (mm/s)', 'FontSize',font_size)




% Correlations
XLF = [v_sv_list(1)*ones(size(PLFA1)); v_sv_list(2)*ones(size(PLFA1)); v_sv_list(3)*ones(size(PLFA1)); ...
    v_sv_list(4)*ones(size(PLFA1)); v_sv_list(5)*ones(size(PLFA1))];
YLF = [PLFA1; PLFA2; PLFA3; PLFA4; PLFA5];
[r_LF, ~] = corr(XLF, YLF)



