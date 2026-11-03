function run_selected(refine)
if nargin<1,refine=false;end
root=paper_data_root();P=design_parameters();
P.outputDt=.1;P.denseDt=.01;P.storeSolutions=true;P.captureMesh=true;P.maxWallSeconds=900;
P.maxStep=.25;P.RelTol=2e-9;P.AbsTol=1e-12;
if refine,P.maxStep=.125;P.RelTol=2e-10;P.AbsTol=[1e-13*ones(60,1);1e-12*ones(P.dim-60,1)];name='selected_refined';else,name='selected';end
folder=fullfile(root,'results',name);if ~exist(folder,'dir'),mkdir(folder);end
C=theorem_certificates(P);write_json(fullfile(folder,'certificate.json'),C);assert(C.analyticLocalPass);
S=audit_core('simulate',P,true);save(fullfile(folder,'trajectory.mat'),'P','S','C','-v7');
audit_final(name,'full');fprintf('SELECTED_DENSE_DONE %s\n',name);
end
