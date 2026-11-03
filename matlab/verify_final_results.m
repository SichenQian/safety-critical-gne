function verify_final_results()
root=paper_data_root();verify_scaling_and_ablation;
report=struct();
for name={'selected','selected_refined','untightened','untightened_refined'}
 n=name{1};folder=fullfile(root,'results',n);R=jsondecode(fileread(fullfile(folder,'final_full.json')));G=jsondecode(fileread(fullfile(folder,'graph_statistics.json')));E=jsondecode(fileread(fullfile(folder,'geometry_extrema.json')));B=jsondecode(fileread(fullfile(folder,'region_check.json')));
 assert(R.validQP==R.samples && max(R.maxKKT)<1e-8 && G.maxComponents==1 && G.minimumBottleneckTreeWeight>0 && B.trajectoryWithinRegion);
 if startsWith(n,'selected')
  D=load(fullfile(folder,'trajectory.mat'),'P');P=D.P;K=jsondecode(fileread(fullfile(folder,'refined_peaks.json')));
  assert(all(R.minimumBarriers>0) && R.minimumPsi>0 && all(R.maxECBF<0) && max([R.maxTracking;K.tracking(:,2)])<P.trackingBudget && max([R.maxDeltaBall;K.deltaBall(:,2)])<P.connMismatchBudget);
  assert(R.minimumLambda2>P.lambdaLower && R.maxInputNorm<P.auditInputRadius && R.maxOptimizerInputNorm<P.auditInputRadius);
  assert(size(R.crossings,1)==2 && isequal(R.crossings(:,2:3),[2 5;3 4]) && isequal(R.crossings(:,7),[1;0]));validate_switches(n);
 else,assert(E.collision.h<0 && E.obstacle.h<0 && E.tree.h<0 && R.minimumLambda2>.04);end
end
for name={'local','local_refined'}
 n=name{1};validate_local_certificate(n);D=load(fullfile(root,'results',n,'trajectory.mat'),'P','S','C');rows=reshape([2*D.C.omittedPairs(:)-1,2*D.C.omittedPairs(:)].',[],1);assert(all(D.S.z(:,D.P.iL(rows))==0,'all'));
 report.(n)=struct('omittedRows',rows,'omittedMultipliersZero',true);
end
for kind={'selected','untightened'}
 n=kind{1};A=load(fullfile(root,'results',n,'trajectory.mat'),'S');B=load(fullfile(root,'results',[n '_refined'],'trajectory.mat'),'S');
 assert(isequal(A.S.t,B.S.t));V=struct('maximumSlowStateDifference',max(abs(A.S.z(:,1:60)-B.S.z(:,1:60)),[],'all'),...
 'maximumFullStateDifference',max(abs(A.S.z-B.S.z),[],'all'),'switchTimeDifference',abs(A.S.sensingCrossings(:,1)-B.S.sensingCrossings(:,1)));
 E1=jsondecode(fileread(fullfile(root,'results',n,'geometry_extrema.json')));E2=jsondecode(fileread(fullfile(root,'results',[n '_refined'],'geometry_extrema.json')));
 V.barrierMinimumDifference=abs([E1.collision.h E1.obstacle.h E1.tree.h]-[E2.collision.h E2.obstacle.h E2.tree.h]);report.(n)=V;
end
write_json(fullfile(root,'tables','final_validation.json'),report);fprintf('FINAL_RESULTS_VALIDATION_PASS\n');
end
