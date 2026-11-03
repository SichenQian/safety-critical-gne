function report=paper_verify_archive(repo,D)
% Separate arithmetic/initialization checks from recorded dense diagnostics.
P=design_parameters();C=theorem_certificates(P);assert(C.theorem1&&C.analyticLocalPass);
for name={'dc','dsafe','lambdaLower','c','gamma','epsilon','xi','K1','K2','DeltaInt','DeltaObs','DeltaC','trackingBudget','connMismatchBudget'}
 assert(isequal(P.(name{1}),D.P.(name{1})),'Parameter/cache mismatch: %s',name{1});
end
assert(isequal(P.y0,D.P.y0)&&isequal(P.ystar,D.P.ystar)&&isequal(P.obs,D.P.obs)&&isequal(P.rObs,D.P.rObs));
report=struct('parameterAndLocalInequalitiesPass',true,'theorem1LowerBounds',[C.cLower C.gammaLower],...
 'tildeL',C.tildeL,'localRadius',[C.rTheta1 C.rFast],'epsilonUpper',C.epsUpper);
for name={'selected','selected_refined','untightened','untightened_refined'}
 n=name{1};folder=fullfile(repo,'reference',n);
 R=jsondecode(fileread(fullfile(folder,'final_full.json')));
 G=jsondecode(fileread(fullfile(folder,'graph_statistics.json')));
 assert(R.validQP==R.samples&&max(R.maxKKT)<1e-8&&G.maxComponents==1&&G.minimumBottleneckTreeWeight>0);
 if startsWith(n,'selected')
  K=jsondecode(fileread(fullfile(folder,'refined_peaks.json')));
  assert(all(R.minimumBarriers>0)&&all(R.maxECBF<0));
  assert(max([R.maxTracking;K.tracking(:,2)])<P.trackingBudget);
  assert(max([R.maxDeltaBall;K.deltaBall(:,2)])<P.connMismatchBudget);
 end
 report.(n)=struct('recordedAuditedSamples',R.samples,'maximumKKT',max(R.maxKKT),'maximumComponents',G.maxComponents);
end
assert(isequal(D.events(:,2:3),[2 5;3 4])&&isequal(D.events(:,7),[1;0]));
I=D.local.initial;assert(I.certified&&I.allDualZero&&I.eTheta<C.rTheta1&&I.eFast<C.rFast);
assert(all(D.local.initialState(P.iL)==0));
assert(D.local.slow(1)<C.rTheta1&&D.local.fast(1)<C.rFast);
report.localInitializationPass=true;report.switchingEvents=D.events;
N=jsondecode(fileread(fullfile(repo,'verification','nominal_hb.json')));
assert(N.stateDimension==60&&N.gainCheck.pass&&N.maximumRhsDiscrepancy<1e-11);
assert(N.parameters.c==P.c&&N.parameters.gamma==P.gamma&&N.parameters.dc==P.dc);
assert(isequal(N.parameters.y0,P.y0)&&isequal(N.parameters.v0,P.v0));
assert(N.standard.minimumLambda2>P.lambdaLower&&N.refined.minimumLambda2>P.lambdaLower);
assert(N.precisionComparison.minimumLambda2Difference<1e-7);
report.nominalMinimumLambda2=[N.standard.minimumLambda2,N.refined.minimumLambda2];
A=theorem2_condition_audit();archived=jsondecode(fileread(fullfile(repo,'verification','theorem2_audit.json')));
assert(audit_matches(jsondecode(jsonencode(A)),archived),'Theorem 2 scope report is stale; regenerate it with theorem2_condition_audit.');
report.theorem2DesignAndTrajectoryChecksPass=true;
report.theorem2RegionalUniformCertificate=A.regionalUniformCertificate;
report.theorem2ContinuousTimeProof=A.continuousTimeProof;
report.scope='Recomputed design/local inequalities and recorded dense-trajectory diagnostics; regional uniform hypotheses are stated in docs/numerical_details.md.';
fprintf('PAPER_PARAMETER_ARCHIVE_VALIDATION_PASS\n');
end

function ok=audit_matches(A,B)
% Preserve scope text exactly and allow small floating-point platform changes.
if isnumeric(A)&&isnumeric(B)
 ok=isequal(size(A),size(B))&&all(abs(A-B)<=5e-15+1e-10*max(abs(A),abs(B)),'all');
elseif isstruct(A)&&isstruct(B)
 fields=fieldnames(A);ok=isequal(size(A),size(B))&&isequal(sort(fields),sort(fieldnames(B)));
 for j=1:numel(A)
  for k=1:numel(fields),if ~ok,return;end;ok=audit_matches(A(j).(fields{k}),B(j).(fields{k}));end
 end
elseif iscell(A)&&iscell(B)
 ok=isequal(size(A),size(B));
 for k=1:numel(A),if ~ok,return;end;ok=audit_matches(A{k},B{k});end
elseif ischar(A)&&ischar(B)&&isfinite(str2double(A))&&isfinite(str2double(B))
 ok=audit_matches(str2double(A),str2double(B));
else
 ok=isequaln(A,B);
end
end
