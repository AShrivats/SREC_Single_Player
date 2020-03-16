%Plotting sample paths
%want to plot implied price dynamics of this model
%close all
load('sg_1f5p.mat')
g_npi = gen_opt;
t_npi = trade_opt;
load('sg_1f5p_pi.mat')
g_pi = gen_opt;
t_pi = trade_opt;
eta_npi = 0;
psi_npi = 0;
eta_pi = eta;
psi_pi = psi;

rp = 1;
b_final_pi = zeros(rp, 1);
gen_final_pi = zeros(rp, 1);
trade_final_pi = zeros(rp, 1);
costs_final_pi = zeros(rp, 1);

b_final_simple_npi = zeros(rp, 1);
gen_final_npi = zeros(rp, 1);
trade_final_npi = zeros(rp, 1);
costs_final_npi = zeros(rp, 1);


S0 = S_grid(62); %initial price S_grid 150
b0 = b_grid(101); % assume the firm starts with nothing
% we assume the firm holds their optimal behaviour for the entire time dt

for sim_num = 1:rp
    sim_noise = normrnd(0, sqrt(dt), [1, time_steps*num_per]) * sigma_f;
    e_noise = normrnd(0, nu*sqrt(dt), [1, time_steps*num_per]);

    % to store relevant parameters for active price impacts
    S_path_npi = zeros(time_steps*num_per+1, 1);
    S_path_pi = zeros(time_steps*num_per+1, 1);

    b_path_npi = zeros(time_steps*num_per+1, 1);
    b_path_pi = zeros(time_steps*num_per+1, 1);

    costs_npi = NaN(time_steps*num_per+1, 1);
    costs_pi = NaN(time_steps*num_per+1, 1);

    gen_npi = NaN(time_steps*num_per+1, 1);
    gen_pi = NaN(time_steps*num_per+1, 1);

    trade_npi = NaN(time_steps*num_per+1, 1);
    trade_pi = NaN(time_steps*num_per+1, 1);

    b_path_npi(1) = b0;
    b_path_pi(1) =b0;
    S_path_npi(1) = S0;
    S_path_pi(1) = S0;

    for i = 1:time_steps*num_per
        [X, Y] = meshgrid(b_grid, S_grid);
        
        n = floor((i-1) / time_steps) + 1;
        t = mod(i-1, time_steps)+1;
        % calculating generation and trading for price impacts active
        g_mat_npi = squeeze(g_npi(n, t, :,:));
        g_mat_pi = squeeze(g_pi(n, t, :,:));

        t_mat_npi = squeeze(t_npi(n, t,:,:));
        t_mat_pi = squeeze(t_pi(n, t,:,:));

        gen_npi(i) = interp2(X,Y,g_mat_npi,b_path_npi(i), S_path_npi(i));
        gen_pi(i) = interp2(X,Y,g_mat_pi,b_path_pi(i), S_path_pi(i));

        trade_npi(i) = interp2(X,Y,t_mat_npi,b_path_npi(i), S_path_npi(i));
        trade_pi(i) = interp2(X,Y,t_mat_pi,b_path_pi(i), S_path_pi(i));

        costs_npi(i) = 1 / 2 * zeta * max(0,(gen_npi(i) - h))^2 *dt + trade_npi(i)*S_path_npi(i)*dt + 1 / 2 * gamma * trade_npi(i)^2 *dt;
        costs_pi(i) = 1 / 2 * zeta * max(0,(gen_pi(i) - h))^2 *dt + trade_pi(i)*S_path_pi(i)*dt + 1 / 2 * gamma * trade_pi(i)^2 *dt;


        if mod(i, time_steps) == 0 && n < num_per           
            b_path_npi(i+1) = min(b_max, max(0, (b_path_npi(i) + max(0, gen_npi(i)*dt + e_noise(i)) + trade_npi(i)*dt) - req));
            b_path_pi(i+1) = min(b_max, max(0, (b_path_pi(i) + max(0, gen_pi(i)*dt + e_noise(i)) + trade_pi(i)*dt) - req));
            
            costs_npi(i) = costs_npi(i) + pen * max(0, req - b_path_npi(i) - max(0, gen_npi(i)*dt + e_noise(i)) - trade_npi(i)*dt);
            costs_pi(i) = costs_pi(i) + pen * max(0, req - b_path_pi(i) - max(0, gen_pi(i)*dt + e_noise(i)) - trade_pi(i)*dt);

        else
            b_path_npi(i+1) = min(b_max, max(0, (b_path_npi(i) + max(0, gen_npi(i)*dt + e_noise(i)) + trade_npi(i)*dt)));
            b_path_pi(i+1) = min(b_max, max(0, (b_path_pi(i) + max(0, gen_pi(i)*dt + e_noise(i)) + trade_pi(i)*dt)));

        end
        S_path_npi(i+1) = max(0, min(pen, S_path_npi(i) + mu_f * dt - psi_npi * gen_npi(i) * dt - psi_npi*e_noise(i) + eta_npi * trade_npi(i) * dt + sim_noise(i)));
        S_path_pi(i+1) = max(0, min(pen, S_path_pi(i) + mu_f * dt - psi_pi * gen_pi(i) * dt - psi_pi*e_noise(i) + eta_pi * trade_pi(i) * dt + sim_noise(i)));

    end
    if b_path_npi(end) >= req
        true_costs_npi = sum(costs_npi(1:end-1));
    else
        true_costs_npi = sum(costs_npi(1:end-1)) + pen * (req - b_path_npi(end));
    end
   
    if b_path_pi(end) >= req
        true_costs_pi = sum(costs_pi(1:end-1));
    else
        true_costs_pi = sum(costs_pi(1:end-1)) + pen * (req - b_path_pi(end));
    end
    

    % summary statistics for active price impacts
%     b_final(sim_num) = b_path(end);
%     gen_final(sim_num) = sum(gen(1:end-1)*dt);
%     trade_final(sim_num) = sum(trade(t:end-1)*dt);
%     costs_final(sim_num) = true_costs*-1;
%     alt2_costs_final(sim_num) = alt2_true_costs*-1;

    sim_num
end
%%
f = figure();
np=6;
subplot(np,1,1);
T_grid_mult = linspace(0, T*num_per, 251);
plot(T_grid_mult,b_path_npi.', 'LineWidth',1)
hold on
plot(T_grid_mult,b_path_pi, 'LineWidth',1)

% hold on
% plot(T_grid, b_path_simple)
title('Banked SRECs')

subplot(np,1,2);
plot(T_grid_mult,gen_npi, 'LineWidth',1)
hold on
plot(T_grid_mult,gen_pi, 'LineWidth',1)

% plot(T_grid,gen_simple, 'LineWidth',1)

title('Generation rate')

subplot(np,1,3);
plot(T_grid_mult,trade_npi, 'LineWidth',1)
hold on
plot(T_grid_mult,trade_pi, 'LineWidth',1)

% plot(T_grid, trade_simple, 'LineWidth', 1)
title('Trading rate')

subplot(np,1,4);
plot(T_grid_mult, S_path_npi, 'LineWidth',1)
hold on
plot(T_grid_mult, S_path_pi, 'LineWidth',1)


% plot(T_grid, S_path_simple, 'LineWidth', 1)
title('SREC price')

subplot(np,1,5);
plot(T_grid_mult,costs_npi, 'LineWidth',1)
hold on
plot(T_grid_mult,costs_pi, 'LineWidth',1)

% plot(T_grid, costs_simple, 'LineWidth', 1)
title('Instantaneous incurred costs')

e_noise_plt = NaN(num_per*time_steps + 1, 1);
e_noise_plt(1:num_per*time_steps) = e_noise;
subplot(np,1,6);
plot(T_grid_mult, cumsum(e_noise_plt), 'LineWidth', 1)
hold on
title('Cumulative production noise')
xlabel('Time', 'fontsize', 14)

%save2pdf("path_multi_p_imp", f, 600)

%SAVE2PDF Saves a figure as a properly cropped pdf
%
%   save2pdf(pdfFileName,handle,dpi)
%
%   - pdfFileName: Destination to write the pdf to.
%   - handle:  (optional) Handle of the figure to write to a pdf.  If
%              omitted, the current figure is used.  Note that handles
%              are typically the figure number.
%   - dpi: (optional) Integer value of dots per inch (DPI).  Sets
%          resolution of output pdf.  Note that 150 dpi is the Matlab
%          default and this function's default, but 600 dpi is typical for
%          production-quality.
%
%   Saves figure as a pdf with margins cropped to match the figure size.

%   (c) Gabe Hoffmann, gabe.hoffmann@gmail.com
%   Written 8/30/2007
%   Revised 9/22/2007
%   Revised 1/14/2007

function save2pdf(pdfFileName,handle,dpi)

% Verify correct number of arguments
error(nargchk(0,3,nargin));

% If no handle is provided, use the current figure as default
if nargin<1
    [fileName,pathName] = uiputfile('*.pdf','Save to PDF file:');
    if fileName == 0; return; end
    pdfFileName = [pathName,fileName];
end
if nargin<2
    handle = gcf;
end
if nargin<3
    dpi = 150;
end

% Backup previous settings
prePaperType = get(handle,'PaperType');
prePaperUnits = get(handle,'PaperUnits');
preUnits = get(handle,'Units');
prePaperPosition = get(handle,'PaperPosition');
prePaperSize = get(handle,'PaperSize');

% Make changing paper type possible
set(handle,'PaperType','<custom>');

% Set units to all be the same
set(handle,'PaperUnits','inches');
set(handle,'Units','inches');

% Set the page size and position to match the figure's dimensions
paperPosition = get(handle,'PaperPosition');
position = get(handle,'Position');
set(handle,'PaperPosition',[0,0,position(3:4)]);
set(handle,'PaperSize',position(3:4));

% Save the pdf (this is the same method used by "saveas")
print(handle,'-dpdf',pdfFileName,sprintf('-r%d',dpi),'-cmyk')

% Restore the previous settings
set(handle,'PaperType',prePaperType);
set(handle,'PaperUnits',prePaperUnits);
set(handle,'Units',preUnits);
set(handle,'PaperPosition',prePaperPosition);
set(handle,'PaperSize',prePaperSize);
end
