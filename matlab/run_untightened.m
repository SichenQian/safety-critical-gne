function run_untightened(refine)
if nargin<1,refine=false;end
root=paper_data_root();D=load(fullfile(root,'results','selected','trajectory.mat'),'P','S');P=D.P;
P.DeltaInt=0;P.DeltaObs=0;P.DeltaC=0;
% Use the exact SAME full initial state, including every fast variable.
P.initialState=D.S.z(1,:).';clear D;
name='untightened';if refine,name='untightened_refined';P.RelTol=2e-10;P.AbsTol=[1e-13*ones(60,1);1e-12*ones(P.dim-60,1)];P.maxStep=.125;end
folder=fullfile(root,'results',name);if ~exist(folder,'dir'),mkdir(folder);end
write_json(fullfile(folder,'parameters.json'),P);
S=audit_core('simulate',P,false);save(fullfile(folder,'trajectory.mat'),'P','S','-v7');
audit_final(name,'full');fprintf('UNTIGHTENED_DONE %s\n',name);
end
