%% Density estimation

rng(5);                                  % random seed
n_samples = 200;                         % number of samples
sigma1 = 1;                              % Rayleigh parameter
rayl1 = raylinv(rand(n_samples, 1), 1 ); % sample Rayleigh

x_analytical = 0 : 0.01 : 4;             % x-axis analytixal
pdf_true = x_analytical./ sigma1^2 .* exp(- x_analytical.^2 / (2*sigma1^2) ); % true Rayleigh pdf

% Proposed method
% function call: Envelope_Statistics_tKDE(pixel_values, bw_factor, kernel_type, return_plot)
[est_pdf, x, data_struct1] = envelope_statistics_tKDE(rayl1, 1, "normal", true);

% plot true pdf in left panel
subplot(1,2,1)
plot(x_analytical, pdf_true,'-.','displayname', 'True PDF ~ Rayl.(1)', 'linewidth', 3, 'color', 'k'); hold on; grid on;




%% gCNR estimation

sigma1 = 1;
sigma2 = 2;
rayl1 = raylinv(rand(n_samples, 1), sigma1);
rayl2 = raylinv(rand(n_samples, 1), sigma2);
rayl1 = [rayl1; zeros(200,1)];
rayl2 = [rayl2; zeros(1800,1)];

% Proposed method
% function call: gCNR_tKDE(pixels_i, pixels_o, return_plot, bw_factor, kernel_type)
[gCNR_value, data_struct, f] = gCNR_tKDE(rayl1, rayl2, true, 0.8, []);


% plot
x_analytical = 0 : 0.01 : 10;
pdf_i = x_analytical./ sigma1^2 .* exp(- x_analytical.^2 / (2*sigma1^2) );
pdf_o = x_analytical./ sigma2^2 .* exp(- x_analytical.^2 / (2*sigma2^2) );
true_gCNR = 1-trapz(x_analytical, min([pdf_i; pdf_o]));
plot(x_analytical, pdf_i,'-.','displayname', 'True PDF ~ Rayl.(1)', 'linewidth', 3, 'color', 'k'); hold on; grid on;
plot(x_analytical, pdf_o,'--','displayname', 'True PDF ~ Rayl.(2)', 'linewidth', 3, 'color', 'k'); hold on; grid on;
xlim([0,6])