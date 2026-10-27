function file = export_connectivity_fig(dataDir)
%EXPORT_CONNECTIVITY_FIG Save Fig. 2(d) as editable native MATLAB graphics.
% Reads the same full-precision CSV data used for the delivered EPS.
% No simulation is run and no numerical data are changed.
if nargin < 1 || isempty(dataDir)
    dataDir = fullfile(fileparts(mfilename('fullpath')), 'connectivity_check');
end
T = readtable(fullfile(dataDir,'tightened_connectivity.csv'));
U = readtable(fullfile(dataDir,'untightened_connectivity.csv'));
f = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 3.5 1.72],'Name','Fig. 2(d): Algebraic connectivity', ...
    'NumberTitle','off','Renderer','painters','WindowStyle','normal', ...
    'MenuBar','figure','ToolBar','figure');
c = onCleanup(@()close(f));
ax = axes(f,'Units','normalized','Position',[.14 .24 .815 .71], ...
    'FontName','Times New Roman','FontSize',8,'LineWidth',.5, ...
    'Box','on','TickDir','in','Layer','top','GridColor',[.88 .88 .88], ...
    'GridAlpha',1,'Tag','connectivity_axes');
hold(ax,'on');
it = T.time_s<=400; iu = U.time_s<=400;
ht=plot(ax,T.time_s(it),T.lambda2(it),'-','Color',[1 0 0], ...
    'LineWidth',.85,'DisplayName','Tightened','Tag','tightened');
hu=plot(ax,U.time_s(iu),U.lambda2(iu),'--','Color',[0 0 1], ...
    'LineWidth',.85,'DisplayName','Untightened','Tag','untightened');
yline(ax,.2,'k:','LineWidth',.7,'HandleVisibility','off','Tag','reference_02');
text(ax,396,.26,'Reference: 0.2','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontName','Times New Roman','FontSize',8);
xlim(ax,[0 400]); ylim(ax,[0 4.2]); xticks(ax,0:100:400); yticks(ax,0:4);
grid(ax,'on');
xlabel(ax,'Time $t$ (s)','Interpreter','latex','FontSize',9);
ylabel(ax,'$\lambda_2(L_{\chi}(y))$','Interpreter','latex','FontSize',9);
legend(ax,[ht hu],'Location','northeast','Orientation','horizontal', ...
    'Box','off','FontName','Times New Roman','FontSize',7.5);
f.UserData=struct('panel','Fig. 2(d)','dataSource', ...
    'tightened_connectivity.csv; untightened_connectivity.csv', ...
    'timeRange',[0 400],'reference',.2,'nativeEditableGraphics',true);
file=fullfile(dataDir,'algebraic_connectivity.fig');
set(f,'Visible','on'); drawnow; savefig(f,file); set(f,'Visible','off');
fprintf('Saved native editable figure: %s\n',file);
end
