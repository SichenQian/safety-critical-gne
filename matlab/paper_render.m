function report=paper_render(D,outputFolder)
repo=fileparts(fileparts(mfilename('fullpath')));if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
report=struct();
for c={'trajectory','safety','connectivity','edge','exp'}
 name=c{1};f=openfig(fullfile(repo,'presentation',[name '.fig']),'new','invisible');
 paper_update_figure(f,name,D);
 assert(isempty(findall(f,'Type','image')));f.Visible='on';
 savefig(f,fullfile(outputFolder,[name '.fig']));
 exportgraphics(f,fullfile(outputFolder,[name '.pdf']),'ContentType','vector');
 exportgraphics(f,fullfile(outputFolder,[name '.png']),'Resolution',240);
 f.PaperPositionMode='auto';f.InvertHardcopy='off';
 print(f,fullfile(outputFolder,[name '.eps']),'-depsc','-painters');
 report.(name)=struct('editable',true,'visibleOnOpen',true,'lineObjects',numel(findall(f,'Type','line')));
 close(f);
end
fprintf('PAPER_FIGURES_RENDERED\n');
end
