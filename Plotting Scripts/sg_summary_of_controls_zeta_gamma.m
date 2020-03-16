% Plotting
close all
time_vec = [50 40 30 20 10 1];

%load('sg_1f1p_fine.mat')
%load('sg_seasonality_1p.mat')
load('sg_zeta_min_gamma_min.mat')
%load('sg_zeta_max_gamma_max.mat')
ctrls= figure();
for tt = 1:numel(time_vec)
    V_plt = squeeze(V(time_vec(tt),:,:));
    g_plt = squeeze(gen_opt(time_vec(tt),:,:));
    t_plt = squeeze(trade_opt(time_vec(tt),:,:));

    %S_low = ceil(numel(S_grid)/6);
    S_mid = ceil(numel(S_grid)/2);
    %S_hi = ceil(numel(S_grid)*5/6);
%     subplot(2, 3, 1)
%     set(gca,'FontSize',14)
%     plot(b_grid, g_plt(S_low,:), 'LineWidth', 1)
%     ylim([-50 1050])
%     %legend(strcat('t=', num2str(time_vec(1))),strcat('t=', num2str(time_vec(2))),strcat('t=', num2str(time_vec(3))),strcat('t=', num2str(time_vec(4))),strcat('t=', num2str(time_vec(5))),strcat('t=', num2str(time_vec(6))), 'Location','northeast')
%     hold on
%     title(strcat("S = ", num2str(round(S_grid(S_low)))))

    subplot(2, 1, 1)
    set(gcf,'Position',[100 100 300 500])
    plot(b_grid, g_plt(S_mid,:), 'LineWidth', 1)
    ylim([-50 2050])
    hold on
    title("$\zeta = 0.2, \gamma = 0.2$", "interpreter", "latex", 'fontsize', 16)
    ylabel("Optimal Generation Rate", 'fontsize', 16)
    set(gca,'FontSize',14)


%     subplot(2, 3, 3)
%     set(gca,'FontSize',14)
%     plot(b_grid, g_plt(S_hi,:), 'LineWidth', 1)
%     ylim([-50 1050])
% 
%     hold on
%     title(strcat("S = ", num2str(round(S_grid(S_hi)))))
% 
% 
%     subplot(2, 3, 4)
%     set(gca,'FontSize',14)
%     plot(b_grid, t_plt(S_low,:), 'LineWidth', 1)
%     ylim([-500 500])
%     hold on
    subplot(2, 1, 2)
    set(gca,'FontSize',14)
    plot(b_grid, t_plt(S_mid,:), 'LineWidth', 1)
    ylim([-800 800])
    ylabel("Optimal Trading Rate", 'fontsize', 16)
    hold on
    xlabel("Banked SRECs", 'fontsize', 16)
    set(gca,'FontSize',14)

%     subplot(2, 3, 6)
%     set(gca,'FontSize',14)
%     plot(b_grid, t_plt(S_hi,:), 'LineWidth', 1)
%     ylim([-500 500])
% 
%     hold on

end
save2pdf("zeta02gamma02", ctrls, 600)
% set(gca, 'color', 'none')
%save2pdf("sg_optimal_behaviour", ctrls, 600)
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