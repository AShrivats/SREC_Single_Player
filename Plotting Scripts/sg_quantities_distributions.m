%Plotting sample paths
%want to plot implied price dynamics of this model
%close all
%load('single_period_50_timesteps.mat')
%load('single_period_h_low.mat')
load('sg_h0.mat')

S0 = S_grid(62); %initial price S_grid 150
b0 = b_grid(1); % assume the firm starts with nothing
% we assume the firm holds their optimal behaviour for the entire time dt
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
  

    % summary statistics for active price impacts
    b_final(sim_num) = b_path(end);
    gen_final(sim_num) = sum(gen(1:end-1)*dt);
    trade_final(sim_num) = sum(trade(t:end-1)*dt);
    costs_final(sim_num) = true_costs*-1;

    sim_num
    
end
%banked SRECs
mean(b_final)
std(b_final)
quantile(b_final, 0.25)
quantile(b_final, 0.75)
quantile(b_final, 0.05)
quantile(b_final, 0.95)
% total planned generation
mean(gen_final)
std(gen_final)
quantile(gen_final, 0.25)
quantile(gen_final, 0.75)
quantile(gen_final, 0.05)
quantile(gen_final, 0.95)
% total trading
mean(trade_final)
std(trade_final)
quantile(trade_final, 0.25)
quantile(trade_final, 0.75)
quantile(trade_final, 0.05)
quantile(trade_final, 0.95)
% profit
mean(costs_final)
std(costs_final)
quantile(costs_final, 0.25)
quantile(costs_final, 0.75)
quantile(costs_final, 0.05)
quantile(costs_final, 0.95)

%plot histograms
fig = figure(1);
axes1 = axes('Parent',fig,...
    'Position',[0.20  0.16 0.70  0.77],...
    'FontSize',18);
box(axes1,'on');
hold(axes1,'all');
% total generation
histogram(gen_final)
xlabel("Total Generated SRECs",'fontsize',24)
ylabel("Frequency",'fontsize',24)
title("Histogram of Total Generated SRECs", 'fontsize', 24)

% total trading
fig2 = figure(2);
axes2 = axes('Parent',fig2,...
    'Position',[0.20  0.16 0.70  0.77],...
    'FontSize',18);
box(axes2,'on');
hold(axes2,'all');
histogram(trade_final)
xlabel("Total Purchased/Sold SRECs",'fontsize',24)
ylabel("Frequency",'fontsize',24)
title("Histogram of Total Purchased/Sold SRECs",'fontsize', 24)

% total trading
fig3 = figure(3);
axes3 = axes('Parent',fig3,...
    'Position',[0.20  0.16 0.70  0.77],...
    'FontSize',18);
box(axes3,'on');
hold(axes3,'all');
histogram(costs_final)
xlabel("Profit",'fontsize',24)
ylabel("Frequency",'fontsize',24)
title("Histogram of Profit",'fontsize', 24)

% total trading
fig4 = figure(4);
axes4 = axes('Parent',fig4,...
    'Position',[0.20  0.16 0.70  0.77],...
    'FontSize',18);
box(axes4,'on');
hold(axes4,'all');
histogram(b_final)
xlabel("Total Accumulated SRECs",'fontsize',24)
ylabel("Frequency",'fontsize',24)
title("Histogram of Total Accumulated SRECs",'fontsize', 24)

save2pdf("gen_hist",fig,600)
save2pdf("trade_hist",fig2,600)
save2pdf("profit_hist",fig3,600)
save2pdf("b_hist",fig4,600)


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
