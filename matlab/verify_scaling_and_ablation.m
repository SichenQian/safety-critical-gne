function verify_scaling_and_ablation()
root=paper_data_root();A=load(fullfile(root,'results','selected','trajectory.mat'),'P','S');B=load(fullfile(root,'results','untightened','trajectory.mat'),'P','S');
P=A.P;Q=B.P;fields=setdiff(union(fieldnames(P),fieldnames(Q)),{'initialState'});changed={};
for k=1:numel(fields),f=fields{k};if ~isequaln(P.(f),Q.(f)),changed{end+1}=f;end,end
assert(isequal(sort(changed),sort({'DeltaInt','DeltaObs','DeltaC'})));
assert(isequal(A.S.z(1,:),B.S.z(1,:))&&Q.DeltaInt==0&&Q.DeltaObs==0&&Q.DeltaC==0);
AR=load(fullfile(root,'results','selected_refined','trajectory.mat'),'P','S');BR=load(fullfile(root,'results','untightened_refined','trajectory.mat'),'P','S');
assert(isequal(AR.S.z(1,:),BR.S.z(1,:)));fr=setdiff(union(fieldnames(AR.P),fieldnames(BR.P)),{'initialState','DeltaInt','DeltaObs','DeltaC'});
for j=1:numel(fr),assert(isequaln(AR.P.(fr{j}),BR.P.(fr{j})));end
er=zeros(1,3);
for k=round(linspace(1,numel(A.S.t),41))
 z=A.S.z(k,:).';D=audit_state(P,z,true);r=D.A*z([P.iU P.iZ])+D.b;
 er(1)=max(er(1),D.pairSumAgreement);
 er(2)=max(er(2),max(abs(r(P.rObsRows)-P.barrierScale*(D.Ac(P.cObs,:)*z(P.iU)+D.bc(P.cObs)))));
 er(3)=max(er(3),abs(sum(D.A(P.rConn,11:end)*z(P.iZ))));
end
assert(max(er)<1e-10);
write_json(fullfile(root,'tables','ablation_and_scaling_validation.json'),struct('changedParameterFields',{changed},'fullInitialStateExactlyEqual',true,...
 'refinedFullInitialStateExactlyEqual',true,'refinedNonmarginParametersExactlyEqual',true,...
 'pairAggregateError',er(1),'obstacleScaleError',er(2),'connectivityMismatchCancellationError',er(3),...
 'alpha1',P.alpha1,'alpha2',P.alpha2,'product',P.alpha1*P.alpha2,'sum',P.alpha1+P.alpha2,...
 'activationLHS',2*P.dc*2*P.regionSpeedMax,'activationRHS',P.alpha1*(P.dc^2-P.dsafe^2),...
 'positiveScalingPreservesAggregatedFeasibility',true,'regularizedOptimizerInvarianceClaimed',false));
fprintf('ABLATION_SCALING_VALIDATION_PASS\n');
end
