%Plotting sample paths
%want to plot implied price dynamics of this model
close all
load('sg_pi_review_baseline.mat')

S0 = S_grid(62); %initial price S_grid 150
b0 = b_grid(40); % assume the firm starts with nothing
% we assume the firm holds their optimal behaviour for the entire time dt
rp = 1;
b_final_min = zeros(rp, 1);
gen_final_min = zeros(rp, 1);
trade_final_min = zeros(rp, 1);
costs_final_min = zeros(rp, 1);

b_final_mid = zeros(rp, 1);
gen_final_mid = zeros(rp, 1);
trade_final_mid = zeros(rp, 1);
costs_final_mid = zeros(rp, 1);

b_final_max =  zeros(rp, 1);
gen_final_max=  zeros(rp, 1);
trade_final_max = zeros(rp, 1);
costs_final_max = zeros(rp, 1);
load('sg_pi_etamin_psimin.mat')
V_min = V;
ge_min = gen_opt;
tr_min = trade_opt;
psi_min = psi;
eta_min = eta;
load('sg_pi_etamid_psimid.mat')
V_mid = V;
ge_mid = gen_opt;
tr_mid = trade_opt;
psi_mid = psi;
eta_mid = eta;
load('sg_pi_etamax_psimax.mat')
V_max = V;
ge_max = gen_opt;
tr_max = trade_opt;
psi_max = psi;
eta_max = eta;

diff_costs = zeros(rp, 1);

g_costs_final = zeros(rp,1);

t_costs_final = zeros(rp, 1);

for sim_num = 1:rp
    sim_noise = normrnd(0, sqrt(dt), [1, time_steps]) * sigma_f;
    e_noise = normrnd(0, nu*sqrt(dt), [1, time_steps]);

    % to store relevant parameters for active price impacts
    S_path_min = zeros(time_steps+1, 1);
    S_path_mid = zeros(time_steps+1, 1);
    S_path_max = zeros(time_steps+1, 1);

    b_path_min = zeros(time_steps+1, 1);
    b_path_mid = zeros(time_steps+1, 1);
    b_path_max = zeros(time_steps+1, 1);

    costs_min = NaN(time_steps+1, 1);
    costs_mid = NaN(time_steps+1, 1);
    costs_max = NaN(time_steps+1, 1);

    gen_min = NaN(time_steps+1, 1);
    gen_mid = NaN(time_steps+1, 1);
    gen_max = NaN(time_steps+1, 1);

    trade_min = NaN(time_steps+1, 1);
    trade_mid = NaN(time_steps+1, 1);
    trade_max = NaN(time_steps+1, 1);

    b_path_min(1) = b0;
    b_path_mid(1) = b0;
    b_path_max(1) = b0;

    S_path_min(1) = S0;
    S_path_mid(1) = S0;
    S_path_max(1) = S0;

    for i = 1:time_steps
        [X, Y] = meshgrid(b_grid, S_grid);
        % calculating generation and trading for price impacts active
        g_min = squeeze(ge_min(i, :,:));
        t_min = squeeze(tr_min(i,:,:));
        gen_firm_min = interp2(X,Y,g_min,b_path_min(i), S_path_min(i));
        trade_firm_min = interp2(X,Y,t_min,b_path_min(i), S_path_min(i));
        gen_min(i) = gen_firm_min;
        trade_min(i) = trade_firm_min;

        costs_min(i) = 1 / 2 * zeta * max(0,(gen_min(i) - h))^2 *dt + trade_min(i)*S_path_min(i)*dt + 1 / 2 * gamma * trade_min(i)^2 *dt;
        b_path_min(i+1) = min(max(0, (b_path_min(i) + max(0, gen_min(i)*dt + e_noise(i)) + trade_min(i)*dt)), b_max);
        S_path_min(i+1) = max(0, min(pen, S_path_min(i) + mu_f * dt - psi_min * max(0, gen_min(i)*dt + e_noise(i)) + eta_min * trade_min(i) * dt + sim_noise(i)));

        g_mid = squeeze(ge_mid(i, :,:));
        t_mid = squeeze(tr_mid(i,:,:));
        gen_firm_mid = interp2(X,Y,g_mid,b_path_mid(i), S_path_mid(i));
        trade_firm_mid = interp2(X,Y,t_mid,b_path_mid(i), S_path_mid(i));
        gen_mid(i) = gen_firm_mid;
        trade_mid(i) = trade_firm_mid;

        costs_mid(i) = 1 / 2 * zeta * max(0,(gen_mid(i) - h))^2 *dt + trade_mid(i)*S_path_mid(i)*dt + 1 / 2 * gamma * trade_mid(i)^2 *dt;
        b_path_mid(i+1) = min(max(0, (b_path_mid(i) + max(0, gen_mid(i)*dt + e_noise(i)) + trade_mid(i)*dt)), b_max);
        S_path_mid(i+1) = max(0, min(pen, S_path_mid(i) + mu_f * dt - psi_mid * max(0, gen_mid(i)*dt + e_noise(i)) + eta_mid * trade_mid(i) * dt + sim_noise(i)));

        
        g_max = squeeze(ge_max(i, :,:));
        t_max = squeeze(tr_max(i,:,:));
        gen_firm_max = interp2(X,Y,g_max,b_path_max(i), S_path_max(i));
        trade_firm_max = interp2(X,Y,t_max,b_path_max(i), S_path_max(i));
        gen_max(i) = gen_firm_max;
        trade_max(i) = trade_firm_max;

        costs_max(i) = 1 / 2 * zeta * max(0,(gen_max(i) - h))^2 *dt + trade_max(i)*S_path_max(i)*dt + 1 / 2 * gamma * trade_max(i)^2 *dt;
        b_path_max(i+1) = min(max(0, (b_path_max(i) + max(0, gen_max(i)*dt + e_noise(i)) + trade_max(i)*dt)), b_max);
        S_path_max(i+1) = max(0, min(pen, S_path_max(i) + mu_f * dt - psi_max * max(0, gen_max(i)*dt + e_noise(i)) + eta_max * trade_max(i) * dt + sim_noise(i)));

    end
    if b_path_min(end) >= req
        true_costs_min = sum(costs_min(1:end-1));
    else
        true_costs_min = sum(costs_min(1:end-1)) + pen * (req - b_path_min(end));
    end
    
    if b_path_mid(end) >= req
        true_costs_mid = sum(costs_mid(1:end-1));
    else
        true_costs_mid = sum(costs_mid(1:end-1)) + pen * (req - b_path_mid(end));
    end
  
    if b_path_max(end) >= req
        true_costs_max = sum(costs_max(1:end-1));
    else
        true_costs_max = sum(costs_max(1:end-1)) + pen * (req - b_path_max(end));
    end

    % summary statistics for active price impacts
    b_final_min(sim_num) = b_path_min(end);
    gen_final_min(sim_num) = sum(gen_min(1:end-1)*dt);
    trade_final_min(sim_num) = sum(trade_min(t:end-1)*dt);
    costs_final_min(sim_num) = true_costs_min*-1;
    
    b_final_mid(sim_num) = b_path_min(end);
    gen_final_mid(sim_num) = sum(gen_min(1:end-1)*dt);
    trade_final_mid(sim_num) = sum(trade_min(t:end-1)*dt);
    costs_final_mid(sim_num) = true_costs_min*-1;
    
    b_final_max(sim_num) = b_path_min(end);
    gen_final_max(sim_num) = sum(gen_min(1:end-1)*dt);
    trade_final_max(sim_num) = sum(trade_min(t:end-1)*dt);
    costs_final_max(sim_num) = true_costs_min*-1;

    sim_num
end

f = figure();
np=6;
subplot(np,1,1);
plot(T_grid,b_path_min.' - req*T_grid, 'LineWidth',1)
hold on
plot(T_grid, b_path_mid.' - req*T_grid,'LineWidth',1)
%plot(T_grid, b_path_max,'LineWidth',1)
title('Banked SRECs relative to pro-rated requirement')

subplot(np,1,2);
plot(T_grid,gen_min, 'LineWidth',1)
hold on
plot(T_grid,gen_mid, 'LineWidth',1)
%plot(T_grid,gen_max, 'LineWidth',1)

title('Generation rate')

subplot(np,1,3);
plot(T_grid,trade_min, 'LineWidth',1)
hold on
plot(T_grid,trade_mid, 'LineWidth',1)
%plot(T_grid,trade_max, 'LineWidth',1)

title('Trading rate')

subplot(np,1,4);
plot(T_grid, S_path_min, 'LineWidth',1)
hold on
plot(T_grid, S_path_mid, 'LineWidth', 1)
%plot(T_grid, S_path_max, 'LineWidth', 1)
title('SREC price')

subplot(np,1,5);
plot(T_grid,costs_min, 'LineWidth',1)
hold on
plot(T_grid, costs_mid, 'LineWidth', 1)
%plot(T_grid, costs_max, 'LineWidth', 1)
title('Instantaneous Incurred Costs')

e_noise_plt = NaN(time_steps + 1, 1);
e_noise_plt(1:time_steps) = e_noise;
subplot(np,1,6);
plot(T_grid, cumsum(e_noise_plt), 'LineWidth', 1)
hold on
title('Cumulative production noise')
xlabel('Time', 'fontsize', 14)

save2pdf("sg_path_active_inactive_alt", f, 600);
% subplot(np,1,6);
% plot(1:time_steps, -psi.*gen + eta.*trade)
% hold on
% title('SREC price drift over time')
% % sum(costs(1:51))

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
