% Plotting
%close all
time_vec = [50 40 30 20 10 1];
% V_imp = V;
% g_imp = gen_opt;
% t_imp = trade_opt;

%load('5p_finegrid_extended_2.mat')
load('sg_1f5p_pi.mat')

ctrls= figure();
for tt = 1:numel(time_vec)
    V_plt = squeeze(V(4, time_vec(tt),:,1:401));
    g_plt = squeeze(gen_opt(4, time_vec(tt),:,1:401));
    t_plt = squeeze(trade_opt(4, time_vec(tt),:,1:401));

    S_low = ceil(numel(S_grid)/6);
    S_mid = ceil(numel(S_grid)/2);
    S_hi = ceil(numel(S_grid)*5/6);
    subplot(2, 3, 1)
    set(gca,'FontSize',14)
    b_grid = b_grid(1:401); %HARD CODED FOR NOW
    plot(b_grid, g_plt(S_low,:), 'LineWidth', 1)
    ylim([-50 1050])
    hold on
    title(strcat("S = ", num2str(round(S_grid(S_low)))))
    ylabel("Optimal Generation Rate", 'fontsize', 14)
    subplot(2, 3, 2)
    set(gca,'FontSize',14)
    plot(b_grid, g_plt(S_mid,:), 'LineWidth', 1)
    ylim([-50 1050])

    hold on
    title(strcat("S = ", num2str(round(S_grid(S_mid)))))

    subplot(2, 3, 3)
    set(gca,'FontSize',14)
    plot(b_grid, g_plt(S_hi,:), 'LineWidth', 1)
    ylim([-50 1050])

    hold on
    title(strcat("S = ", num2str(round(S_grid(S_hi)))))


    subplot(2, 3, 4)
    set(gca,'FontSize',14)
    plot(b_grid, t_plt(S_low,:), 'LineWidth', 1)
    ylim([-500 500])
    hold on
    ylabel("Optimal Trading Rate", 'fontsize', 14)
    subplot(2, 3, 5)
    set(gca,'FontSize',14)
    plot(b_grid, t_plt(S_mid,:), 'LineWidth', 1)
    ylim([-500 500])

    hold on
    xlabel("Banked SRECs", 'fontsize', 14)

    subplot(2, 3, 6)
    set(gca,'FontSize',14)
    plot(b_grid, t_plt(S_hi,:), 'LineWidth', 1)
    ylim([-500 500])

    hold on

end
%save2pdf("p5o5", ctrls, 600)
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
save2pdf("p_imp_p4o5", ctrls, 600)

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