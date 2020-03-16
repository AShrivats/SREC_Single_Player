%Plotting sample paths
%want to plot implied price dynamics of this model
%close all
load('sg_1f5p_pi.mat')
%load('single_period_h_low.mat')

rp = 1;
b_final = zeros(rp, 1);
gen_final = zeros(rp, 1);
trade_final = zeros(rp, 1);
costs_final = zeros(rp, 1);

b_final_simple = zeros(rp, 1);
gen_final_simple = zeros(rp, 1);
trade_final_simple = zeros(rp, 1);
costs_final_simple = zeros(rp, 1);


% strategy 2


diff_costs = zeros(rp, 1);

g_costs_final = zeros(rp,1);

t_costs_final = zeros(rp, 1);
S0 = S_grid(62); %initial price S_grid 150
b0 = b_grid(1); % assume the firm starts with nothing
b_vec = [b_grid(1) b_grid(101) b_grid(201)];
% we assume the firm holds their optimal behaviour for the entire time dt

for sim_num = 1:rp
    sim_noise = normrnd(0, sqrt(dt), [1, time_steps*num_per]) * sigma_f;
    e_noise = normrnd(0, nu*sqrt(dt), [1, time_steps*num_per]);

    % to store relevant parameters for active price impacts
    S_path = zeros(time_steps*num_per+1, 1);
    b_path1 = zeros(time_steps*num_per+1, 1);
    b_path2 = zeros(time_steps*num_per+1, 1);
    b_path3 = zeros(time_steps*num_per+1, 1);

    costs1 = NaN(time_steps*num_per+1, 1);
    costs2 = NaN(time_steps*num_per+1, 1);
    costs3 = NaN(time_steps*num_per+1, 1);

    gen1 = NaN(time_steps*num_per+1, 1);
    gen2 = NaN(time_steps*num_per+1, 1);
    gen3 = NaN(time_steps*num_per+1, 1);

    trade1 = NaN(time_steps*num_per+1, 1);
    trade2 = NaN(time_steps*num_per+1, 1);
    trade3 = NaN(time_steps*num_per+1, 1);
    b_path1(1) = b_vec(1);
    b_path2(1) = b_vec(2);
    b_path3(1) = b_vec(3);
    
    S_path(1) = S0;
    for i = 1:time_steps*num_per
        [X, Y] = meshgrid(b_grid, S_grid);
        
        n = floor((i-1) / time_steps) + 1;
        t = mod(i-1, time_steps)+1;
        % calculating generation and trading for price impacts active
        g_mat = squeeze(gen_opt(n, t, :,:));
        t_mat = squeeze(trade_opt(n, t,:,:));
        gen1(i) = interp2(X,Y,g_mat,b_path1(i), S_path(i));
        gen2(i) = interp2(X,Y,g_mat,b_path2(i), S_path(i));
        gen3(i) = interp2(X,Y,g_mat,b_path3(i), S_path(i));

        trade1(i) = interp2(X,Y,t_mat,b_path1(i), S_path(i));
        trade2(i) = interp2(X,Y,t_mat,b_path2(i), S_path(i));
        trade3(i) = interp2(X,Y,t_mat,b_path3(i), S_path(i));

        costs1(i) = 1 / 2 * zeta * max(0, (gen1(i) - h))^2 *dt + trade1(i)*S_path(i)*dt + 1 / 2 * gamma * trade1(i)^2 *dt;
        costs2(i) = 1 / 2 * zeta * max(0, (gen2(i) - h))^2 *dt + trade2(i)*S_path(i)*dt + 1 / 2 * gamma * trade2(i)^2 *dt;
        costs3(i) = 1 / 2 * zeta * max(0, (gen3(i) - h))^2 *dt + trade3(i)*S_path(i)*dt + 1 / 2 * gamma * trade3(i)^2 *dt;

    
        if mod(i, time_steps) == 0 && n < num_per  
            b_path1(i+1) = min(b_max, max(0, (b_path1(i) + max(0, gen1(i)*dt + e_noise(i)) + trade1(i)*dt) - req));
            b_path2(i+1) = min(b_max, max(0, (b_path2(i) + max(0, gen2(i)*dt + e_noise(i)) + trade2(i)*dt) - req));
            b_path3(i+1) = min(b_max, max(0, (b_path3(i) + max(0, gen3(i)*dt + e_noise(i)) + trade3(i)*dt) - req));

            costs1(i) = costs1(i) + pen * max(0, req - b_path1(i) - max(0, gen1(i)*dt + e_noise(i)) - trade1(i)*dt);
            costs2(i) = costs2(i) + pen * max(0, req - b_path2(i) - max(0, gen2(i)*dt + e_noise(i)) - trade2(i)*dt);
            costs3(i) = costs3(i) + pen * max(0, req - b_path3(i) - max(0, gen3(i)*dt + e_noise(i)) - trade3(i)*dt);

        else
            b_path1(i+1) = min(b_max, max(0, (b_path1(i) + max(0, gen1(i)*dt + e_noise(i)) + trade1(i)*dt)));
            b_path2(i+1) = min(b_max, max(0, (b_path2(i) + max(0, gen2(i)*dt + e_noise(i)) + trade2(i)*dt)));
            b_path3(i+1) = min(b_max, max(0, (b_path3(i) + max(0, gen3(i)*dt + e_noise(i)) + trade3(i)*dt)));

        end
        S_path(i+1) = max(0, min(pen, S_path(i) + mu_f * dt + sim_noise(i))); % as eta and psi are 0
    %e_noise(45) = -11;
    end
    if b_path1(end) >= req
        true_costs = sum(costs1(1:end-1));
    else
        true_costs = sum(costs1(1:end-1)) + pen * (req - b_path1(end));
    end
   
    if b_path2(end) >= req
        true_costs2 = sum(costs2(1:end-1));
    else
        true_costs2 = sum(costs2(1:end-1)) + pen * (req - b_path2(end));
    end
    
    if b_path3(end) >= req
        true_costs3 = sum(costs3(1:end-1));
    else
        true_costs3 = sum(costs3(1:end-1)) + pen * (req - b_path3(end));
    end

    % summary statistics for active price impacts
%     b_final(sim_num) = b_path(end);
%     gen_final(sim_num) = sum(gen(1:end-1)*dt);
%     trade_final(sim_num) = sum(trade(t:end-1)*dt);
%     costs_final(sim_num) = true_costs*-1;
%     alt2_costs_final(sim_num) = alt2_true_costs*-1;

    sim_num
end
%mean(b_final)


% figure()
% plot(gen_final, trade_final, 'o')
% title("Total Generation vs Total Trading")
% xlabel("Total Generation")
% ylabel("Total Trading")

%comparison between Mean Behaviour and Output
% figure()
% histogram(costs_final, 'facealpha', 0.5)
% hold on
% histogram(alt2_costs_final, 'facealpha', 0.5)
% legend('Output Strategy', 'Mean Behaviour Strategy')
% xlabel("Total Profit")
% ylabel("Frequency")
% title("Comparison of Strategies")
% 
% figure()
% histogram(costs_final - alt2_costs_final)
% xlabel("Difference in Profit (Output Strategy - Mean Behaviour Strategy")
% ylabel("Frequency")
% title("Comparison of Strategies")
%plot histograms
% f1 = figure();
% % total generation
% histogram(gen_final)
% xlabel("Total Generated SRECs")
% ylabel("Frequency")
% title("Histogram of Total Generated SRECs (1000 simulations)")
% 
% % total trading
% f2 = figure();
% histogram(trade_final)
% xlabel("Total Purchased/Sold SRECs")
% ylabel("Frequency")
% title("Histogram of Total Purchased/Sold SRECs (1000 simulations)")
% 
% % costs
% f3 = figure();
% histogram(costs_final)
% xlabel("Total Profit")
% ylabel("Frequency")
% title("Histogram of Total Profit (1000 simulations)")
% yL = get(gca,'YLim');
% l1 = line([mean(costs_final) mean(costs_final)],yL,'Color','r');% np = 5;
% l2 = line([V(1, 62, 1)*-1 V(1,62,1)*-1],yL,'Color','b');% np = 5;
% legend([l1 l2],'Mean of Simulations','Theoretical Mean')

% % figure()
% % histogram(alt1_costs_final)
% % xlabel("Total Revenue")
% % ylabel("Frequency")
% % title("Histogram of Total Revenue (1000 simulations)")
% % mean(costs_final)
% % mean(alt1_costs_final)
% % figure()
% % hist(diff_costs, 20)
% % xlabel("Difference Between Optimal and Alternative Strategies")
% % ylabel("Frequency")
% % title("Histogram of Cost Differences Between Strategies")

f = figure();
np=6;
subplot(np,1,1);
T_grid_mult = linspace(0, T*num_per, 251);
plot(T_grid_mult,b_path1, 'LineWidth',1)
hold on
plot(T_grid_mult,b_path2, 'LineWidth',1)
plot(T_grid_mult,b_path3, 'LineWidth',1)

% hold on
% plot(T_grid, b_path_simple)
title('Banked SRECs')

subplot(np,1,2);
plot(T_grid_mult,gen1, 'LineWidth',1)
hold on
plot(T_grid_mult,gen2, 'LineWidth',1)
plot(T_grid_mult,gen3, 'LineWidth',1)

% plot(T_grid,gen_simple, 'LineWidth',1)

title('Generation rate')

subplot(np,1,3);
plot(T_grid_mult,trade1, 'LineWidth',1)
hold on
plot(T_grid_mult,trade2, 'LineWidth',1)
plot(T_grid_mult,trade3, 'LineWidth',1)

% plot(T_grid, trade_simple, 'LineWidth', 1)
title('Trading rate')

subplot(np,1,4);
plot(T_grid_mult, S_path, 'LineWidth',1)
% hold on
% plot(T_grid, S_path_simple, 'LineWidth', 1)
title('SREC price')

subplot(np,1,5);
plot(T_grid_mult,costs1, 'LineWidth',1)
hold on
plot(T_grid_mult,costs2, 'LineWidth',1)
plot(T_grid_mult,costs3, 'LineWidth',1)

% plot(T_grid, costs_simple, 'LineWidth', 1)
title('Instantaneous incurred costs')

%save2pdf("path_active_inactive", f, 600);
e_noise_plt = NaN(num_per*time_steps + 1, 1);
e_noise_plt(1:num_per*time_steps) = e_noise;
subplot(np,1,6);
plot(T_grid_mult, cumsum(e_noise_plt), 'LineWidth', 1)
hold on
title('Cumulative production noise')
xlabel('Time', 'fontsize', 14)

%save2pdf("5_per_path", f, 600)


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
