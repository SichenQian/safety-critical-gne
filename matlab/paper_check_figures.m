function report=paper_check_figures(folder,D)
% Check numerical correspondence by updating an in-memory copy and comparing.
report=struct();
for c={'trajectory','safety','connectivity','edge','exp'}
 name=c{1};f=openfig(fullfile(folder,[name '.fig']),'new');
 assert(strcmp(f.Visible,'on'));f.Visible='off';
 if strcmp(name,'trajectory')&&(~isstruct(f.UserData)||~isfield(f.UserData,'paperBinding'))
  paper_bind_template(f,name,D);
 end
 before=signature(f);paper_update_figure(f,name,D);after=signature(f);
 assert(numel(before)==numel(after));worst=0;
 for k=1:numel(before)
  assert(isequal(size(before{k}),size(after{k})));
  worst=max(worst,max(abs(before{k}-after{k}),[],'all'));
 end
 assert(worst<1e-10,'Figure %s does not match the saved numerical data (%.3g).',name,worst);
 assert(isempty(findall(f,'Type','image')));
 report.(name)=struct('maximumDataDifference',worst,'editableLineCount',numel(findall(f,'Type','line')),'visibleOnOpen',true);
 close(f);
end
fprintf('PAPER_FIGURE_DATA_VALIDATION_PASS\n');
end
function s=signature(f)
h=findall(f,'Type','line');s=cell(numel(h),1);
for k=1:numel(h),s{k}=[h(k).XData(:) h(k).YData(:)];end
for h=findall(f,'Type','patch').',s{end+1}=[h.XData(:) h.YData(:)];end
for h=findall(f,'Type','rectangle').',s{end+1}=h.Position;end
end
