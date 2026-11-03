function result=reproduce_paper(mode,outputFolder)
%REPRODUCE_PAPER Reproduce the current manuscript's five numerical figures.
% reproduce_paper                  Redraw from the included compact data.
% reproduce_paper('full')          Rerun nominal HB and six filtered/local runs.
% reproduce_paper('verify')        Check the included data and figure archive.
% reproduce_paper('download',dir)  Download all full reference MAT/CSV data.
% reproduce_paper('verify',dir)    Check downloaded or freshly simulated data.
% reproduce_paper('nominal',dir)   Reproduce Table I's nominal HB connectivity.
if nargin<1||isempty(mode),mode='figures';end
assert(ismember(mode,{'figures','full','verify','download','nominal'}),'Unknown reproduction mode.');
repo=fileparts(mfilename('fullpath'));addpath(fullfile(repo,'matlab'));
if nargin<2||isempty(outputFolder)
 outputFolder=fullfile(repo,'runs',[mode '_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
if strcmp(mode,'download'),result=download_reference_data(outputFolder);return;end
if strcmp(mode,'nominal')
 result=run_nominal_hb(outputFolder);
 fprintf('REPRODUCTION_PASS mode=nominal\nOutput: %s\n',outputFolder);return;
end
oldRoot=paper_data_root();restore=onCleanup(@()paper_data_root(oldRoot)); %#ok<NASGU>
for d={'results','tables','figures','verification','manuscript','data'}
 if ~exist(fullfile(outputFolder,d{1}),'dir'),mkdir(fullfile(outputFolder,d{1}));end
end
if strcmp(mode,'full')
 assert(license('test','Optimization_Toolbox'),'Full simulation requires Optimization Toolbox.');
 paper_data_root(outputFolder);
 run_nominal_hb(outputFolder);
 validate_generalization;validate_adopted_equations;validate_error_coordinates;
 run_selected(false);run_selected(true);run_untightened(false);run_untightened(true);
 run_local('local',false);run_local('local_refined',true);
 for name={'selected','selected_refined','untightened','untightened_refined'},final_diagnostics(name{1});end
 D=paper_prepare_data(outputFolder);
 fits=struct('slowWindow',[20 120],'slowSlope',D.slowFit(1),'fastWindow',[0 5*D.P.epsilon/D.P.xi],...
  'fastSlope',D.fastFit(1),'theoreticalNormRateLower',D.local.certificate.rateLower,...
  'fastMechanism','Numerical slope of the initial layer; not a region-uniform rate.');
 write_json(fullfile(outputFolder,'tables','local_fits.json'),fits);
 build_final_tables;verify_final_results;
 save(fullfile(outputFolder,'data','paper_plot_data.mat'),'D','-v7');
elseif strcmp(mode,'verify')&&exist(fullfile(outputFolder,'results','selected','trajectory.mat'),'file')
 paper_data_root(outputFolder);verify_final_results;D=paper_prepare_data(outputFolder);
else
 cache=load(fullfile(repo,'data','paper_plot_data.mat'),'D');D=cache.D;
end
result=paper_verify_archive(repo,D);
if strcmp(mode,'verify')
 result.figures=paper_check_figures(fullfile(repo,'figures'),D);
else
 result.exports=paper_render(D,fullfile(outputFolder,'figures'));
 result.figures=paper_check_figures(fullfile(outputFolder,'figures'),D);
end
result.mode=mode;result.outputFolder=outputFolder;result.matlabVersion=version;
write_json(fullfile(outputFolder,'verification','reproduction_summary.json'),result);
fprintf('REPRODUCTION_PASS mode=%s\nOutput: %s\n',mode,outputFolder);
end
