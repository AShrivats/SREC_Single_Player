%Plotting sample paths
%want to plot implied price dynamics of this model
close all
%load('single_period_50_timesteps.mat')
%load('single_period_h_low.mat')
%load('sg_test.mat')
load('sg_pi_etamin_psimin.mat')

S0 = S_grid(62); %initial price S_grid 150
b0 = b_grid(1); % assume the firm starts with nothing
% we assume the firm holds their optimal behaviour for the entire time dt
rp = 1;
rp = 1000;

b_final = zeros(rp, 1);
gen_final = zeros(rp, 1);
trade_final = zeros(rp, 1);
costs_final = zeros(rp, 1);

b_final_simple = zeros(rp, 1);
gen_final_simple = zeros(rp, 1);
trade_final_simple = zeros(rp, 1);
costs_final_simple = zeros(rp, 1);


% strategy 2
c_gen = 605;
c_trd = -145;
alt1_costs_final = zeros(rp, 1);
alt2_costs_final = zeros(rp, 1);
alt3_costs_final = zeros(rp, 1);
alt4_costs_final = zeros(rp, 1);

diff_costs = zeros(rp, 1);

g_costs_final = zeros(rp,1);

t_costs_final = zeros(rp, 1);

for sim_num = 1:rp
    sim_noise = normrnd(0, sqrt(dt), [1, time_steps]) * sigma_f;
    e_noise = normrnd(0, nu*sqrt(dt), [1, time_steps]);
    % to store relevant parameters for active price impacts
    S_path = zeros(time_steps+1, 1);
    b_path = zeros(time_steps+1, 1);
    costs = NaN(time_steps+1, 1);
    gen = NaN(time_steps+1, 1);
    trade = NaN(time_steps+1, 1);
    
    alt2_costs = NaN(time_steps+1, 1);
    alt2_b_path = zeros(time_steps+1, 1);
    
    b_path(1) = b0;
    S_path(1) = S0;
    for i = 1:time_steps
        [X, Y] = meshgrid(b_grid, S_grid);
        % calculating generation and trading for price impacts active
        g_mat = squeeze(gen_opt(i, :,:));
        t_mat = squeeze(trade_opt(i,:,:));
        gen_firm = interp2(X,Y,g_mat,b_path(i), S_path(i));
        trade_firm = interp2(X,Y,t_mat,b_path(i), S_path(i));
        gen(i) = gen_firm;
        trade(i) = trade_firm;

        costs(i) = 1 / 2 * zeta * (max(0, gen(i) - h))^2 *dt + trade(i)*S_path(i)*dt + 1 / 2 * gamma * trade(i)^2 *dt;

        b_path(i+1) = min(max(0, (b_path(i) + max(0, gen(i)*dt + e_noise(i)) + trade(i)*dt)), b_max);
        S_path(i+1) = max(0, min(pen, S_path(i) + mu_f * dt - psi * max(0, gen(i)*dt + e_noise(i)) + eta * trade(i) * dt + sim_noise(i)));
    end
    if b_path(end) >= req
        true_costs = sum(costs(1:end-1));
    else
        true_costs = sum(costs(1:end-1)) + pen * (req - b_path(end));
    end
   
    if alt2_b_path(end) >= req
        alt2_true_costs = sum(alt2_costs(1:end-1));
    else
        alt2_true_costs = sum(alt2_costs(1:end-1)) + pen * (req - alt2_b_path(end));
    end

    % summary statistics for active price impacts
    b_final(sim_num) = b_path(end);
    gen_final(sim_num) = sum(gen(1:end-1)*dt);
    trade_final(sim_num) = sum(trade(t:end-1)*dt);
    costs_final(sim_num) = true_costs*-1;
    alt2_costs_final(sim_num) = alt2_true_costs*-1;

    diff_costs(sim_num) = costs_final(sim_num) - alt1_costs_final(sim_num);
    sim_num
    
end
mb = mean(b_final)
sb = std(b_final)
q1b = quantile(b_final, 0.25)
q3b = quantile(b_final, 0.75)
skb = skewness(b_final)
kb = kurtosis(b_final)

mg = mean(gen_final)
sg = std(gen_final)
q1g = quantile(gen_final, 0.25)
q3g = quantile(gen_final, 0.75)
skg = skewness(gen_final)
kg = kurtosis(gen_final)

mt = mean(trade_final)
st = std(trade_final)
q1t = quantile(trade_final, 0.25)
q3t = quantile(trade_final, 0.75)
skt = skewness(trade_final)
kt = kurtosis(trade_final)

mc = mean(costs_final)
st = std(costs_final)
q1c = quantile(costs_final, 0.25)
q3c = quantile(costs_final, 0.75)
skc = skewness(costs_final)
kc = kurtosis(costs_final)
%mean(b_final)
% 
% %f = figure('Position', [10 10 600 900]);
% f = figure();
% set(gcf,'Position',[100 100 750 750])
% np=6;
% subplot(np,1,1);
% plot(T_grid,b_path, 'LineWidth',1)
% set(gca,'FontSize',14)
% % hold on
% % plot(T_grid, b_path_simple)
% title('Banked SRECs')
% 
% subplot(np,1,2);
% plot(T_grid,gen, 'LineWidth',1)
% set(gca,'FontSize',14)
% % hold on
% % plot(T_grid,gen_simple, 'LineWidth',1)
% 
% title('Generation rate')
% 
% subplot(np,1,3);
% plot(T_grid,trade, 'LineWidth',1)
% set(gca,'FontSize',14)
% % hold on
% % plot(T_grid, trade_simple, 'LineWidth', 1)
% title('Trading rate')
% 
% subplot(np,1,4);
% plot(T_grid, S_path, 'LineWidth',1)
% set(gca,'FontSize',14)
% % hold on
% % plot(T_grid, S_path_simple, 'LineWidth', 1)
% title('SREC price')
% 
% subplot(np,1,5);
% plot(T_grid,costs, 'LineWidth',1)
% set(gca,'FontSize',14)
% % hold on
% % plot(T_grid, costs_simple, 'LineWidth', 1)
% title('Instantaneous incurred costs')
% 
% e_noise_plt = NaN(time_steps+1, 1);
% e_noise_plt(1:time_steps) = cumsum(e_noise);
% subplot(np,1,6);
% plot(T_grid,e_noise_plt, 'LineWidth',1)
% set(gca,'FontSize',14)
% % hold on
% % plot(T_grid, costs_simple, 'LineWidth', 1)
% title('Cumulative production noise')
% xlabel('Time', 'fontsize', 20)

%save2pdf("sg_path", f, 600);
% subplot(np,1,6);
% plot(1:time_steps, -psi.*gen + eta.*trade)
% hold on
% title('SREC price drift over time')
% % sum(costs(1:51))

% covariation plots
% f1 = figure();
% 
% %S and g
% S_diff = S_path(2:end) - S_path(1:end-1);
% S_diff = S_diff(1:end-1);
% g_diff = gen(2:end) - gen(1:end-1);
% g_diff = g_diff(1:end-1);
% t_diff = trade(2:end) - trade(1:end-1);
% t_diff = t_diff(1:end-1);
% cov_sg = S_diff .* g_diff;
% cov_gt = g_diff .* t_diff;
% cov_st = S_diff .* t_diff;
% subplot(3, 1, 1)
% plot(1:time_steps-1, cumsum(cov_sg))
% title("Covariation between S and Planned Generation")
% subplot(3, 1, 2)
% plot(1:time_steps-1, cumsum(cov_gt))
% title("Covariation between Trading and Planned Generation")
% subplot(3, 1, 3)
% plot(1:time_steps-1, cumsum(cov_st))
% title("Covariation between S and Trading")
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
