function [gCNR_value, data_struct, f] = gCNR_tKDE(pixels_i, pixels_o, return_plot, bw_factor, kernel_type)
    
%% Calculate gCNR from two pixel populations (region "i" and "o") using the transform KDE method.
% 

    arguments
        pixels_i
        pixels_o
        return_plot = [];
        bw_factor = [];
        kernel_type = [];
    end
    
    % in case arguments are left empty
    if isempty(return_plot); return_plot = false; end;
    if isempty(bw_factor);  bw_factor = 0.7;          end;
    if isempty(kernel_type); kernel_type = "normal";        end; 
    
    %% Check and initialize
    data_struct = {};

    % check NaNs
    if anynan(pixels_i); error('NaN in pixels_i'); end
    if anynan(pixels_o); error('NaN in pixels_o'); end

    % make sure values are not complex
    pixels_i = abs(pixels_i(:));
    pixels_o = abs(pixels_o(:));

    %% Handle true zeros
    % We consider true zeros to be special (like delta distributions), 
    % so we take them out and compare them directly later.
    N_i = length(pixels_i);
    N_zeros_i = sum(pixels_i == 0);
    pixels_i = pixels_i(pixels_i ~= 0);

    N_o = length(pixels_o);
    N_zeros_o = sum(pixels_o == 0);
    pixels_o = pixels_o(pixels_o ~= 0);

    %% Estiamate gCNR
    
    % Estimate PDFs (via data in structs)
    [~, x_i, data_struct_i] = Envelope_Statistics_tKDE(pixels_i, bw_factor, kernel_type, false);
    [~, x_o, data_struct_o] = Envelope_Statistics_tKDE(pixels_o, bw_factor, kernel_type, false);

    % Common domain sampling
    x = unique([ x_i; x_o ]);

    pdf_i = backtransform(data_struct_i.x_BC, data_struct_i.pdf_BC, x, data_struct_i.alpha_BC);
    pdf_o = backtransform(data_struct_o.x_BC, data_struct_o.pdf_BC, x, data_struct_o.alpha_BC);

    % Normalize PDFs, might be needed if some density is assigned very close to zero
    pdf_i = pdf_i ./ trapz(x, pdf_i);
    pdf_o = pdf_o ./ trapz(x, pdf_o);

    % Integration of OVL
    
    OVL = trapz(x, min( [pdf_i*(1-N_zeros_i/N_i), ...     % pdf excl. zeros
                         pdf_o*(1-N_zeros_o/N_o)].' )) ...% pdf excl. zeros
          + min(N_zeros_i/N_i, N_zeros_o/N_o);            % compare zeros separately

    gCNR_value = 1 - OVL;

    %% Put in struct
    data_struct.str_i = data_struct_i;
    data_struct.str_o = data_struct_o;
    data_struct.gCNR = gCNR_value; % transformed domain
    data_struct.x = x;
    data_struct.pdf_i = pdf_i;
    data_struct.pdf_o = pdf_o;

    %% Plot if 'return_plot' is true
    if return_plot == true
        f = figure();
    
        hhi = area(x, pdf_i, 'LineStyle','none', 'displayname', '\it {f}_i \rm estimate'); hold on; grid on;
        hho = area(x, pdf_o, 'LineStyle','none', 'displayname', '\it {f}_o\rm estimate');
        hh = area(x, min([pdf_i, pdf_o].'), 'LineStyle','none', 'displayname', 'OVL');
        hhi.FaceColor = 'r';      
        hhi.FaceAlpha = 0.8;      
        hho.FaceColor = 'b';
        hho.FaceAlpha = 0.8;
        hh.FaceColor = [0.6 0.6 0.6];      
        ylabel('Probability density');
        
        xlabel('\itx');
        
        legend();
        plot(NaN,NaN,'display',strcat( "Zeros in f_i: ", num2str(100*N_zeros_i/N_i, 3) ,"%"), 'linestyle', 'none')
        plot(NaN,NaN,'display',strcat( "Zeros zeros in f_o: ", num2str(100*N_zeros_o/N_o, 3) ,"%"), 'linestyle', 'none')
    else
        f=0;
    end

end




