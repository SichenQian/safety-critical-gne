function report=run_nominal_hb(outputFolder)
%RUN_NOMINAL_HB Independently reproduce Table I's nominal-flow lambda2.
% No safety-filter or fast states are integrated. No figures are produced.
if nargin<1||isempty(outputFolder)
 repo=fileparts(fileparts(mfilename('fullpath')));
 outputFolder=fullfile(repo,'runs',['nominal_' datestr(now,'yyyymmdd_HHMMSS')]);
end
P=design_parameters();C=theorem_certificates(P);assert(C.theorem1);
X0=repmat(reshape(P.y0.',1,10),5,1)+P.estimateOffset;
z0=[reshape(X0.',50,1);reshape(P.v0.',10,1)];
zs=[reshape(repmat(reshape(P.ystar.',1,10),5,1).',50,1);zeros(10,1)];
% Check the independently written RHS against the shared model equations.
probes=[z0,zs+1e-3*sin((1:60)'),z0+1e-3*cos((1:60)')];rhsDifference=0;
for k=1:size(probes,2)
 z=probes(:,k);X=reshape(z(1:50),10,5).';V=reshape(z(51:60),2,5).';
 [~,~,un,aux]=audit_core('rows',P,X,V,false);
 expected=[reshape(aux.Xdot.',50,1);reshape(aux.f0.',10,1)+un];
 rhsDifference=max(rhsDifference,norm(nominal_hb_rhs(P,z)-expected,inf));
 assert(norm(nominal_hb_rhs(P,z)-expected,inf)<1e-11*max(1,norm(expected,inf)));
end
[standard,A]=solve_case(P,z0,zs,outputFolder,false);
[refined,B]=solve_case(P,z0,zs,outputFolder,true);
assert(isequal(A.t,B.t));
comparison=struct('maximumSlowStateDifference',max(abs(A.z-B.z),[],'all'),...
 'maximumOutputLambda2Difference',max(abs(A.lambda2-B.lambda2)),...
 'minimumLambda2Difference',abs(standard.minimumLambda2-refined.minimumLambda2));
assert(comparison.minimumLambda2Difference<1e-7);
report=struct('equations','nominal_HB_stacked','stateDimension',60,...
 'stateOrder','col{reshape(X transpose,50,1),reshape(V transpose,10,1)}',...
 'parameters',struct('c',P.c,'gamma',P.gamma,'dc',P.dc,'lambdaLower',P.lambdaLower,...
 'B',P.B,'p',P.p,'y0',P.y0,'v0',P.v0,'initialEstimates',X0,'target',P.ystar),...
 'gainCheck',struct('tildeL',C.tildeL,'cLower',C.cLower,'gammaLower',C.gammaLower,'pass',C.theorem1),...
 'maximumRhsDiscrepancy',rhsDifference,'standard',standard,'refined',refined,...
 'precisionComparison',comparison,'matlabVersion',version,...
 'scope','Finite-horizon numerical trajectory and interpolant diagnostics; not an analytic continuous-time lower bound.');
folder=fullfile(outputFolder,'verification');if ~exist(folder,'dir'),mkdir(folder);end
write_json(fullfile(folder,'nominal_hb.json'),report);
fprintf('NOMINAL_HB_VALIDATION_PASS minimum_lambda2=%.15g t=%.9g\n',standard.minimumLambda2,standard.minimumTime);
end

function [R,S]=solve_case(P,z0,zs,root,refine)
P.RelTol=2e-9;P.AbsTol=1e-12;P.maxStep=.25;P.outputDt=.1;P.denseDt=.01;
name='nominal';if refine,P.RelTol=2e-10;P.AbsTol=1e-13;P.maxStep=.125;name='nominal_refined';end
folder=fullfile(root,'results',name);if ~exist(folder,'dir'),mkdir(folder);end
opts=odeset('RelTol',P.RelTol,'AbsTol',P.AbsTol,'MaxStep',P.maxStep,...
 'Events',@(t,z)lower_event(P,z));
sol=ode15s(@(t,z)nominal_hb_rhs(P,z),[0 P.T],z0,opts);
assert(sol.x(end)==P.T,'Nominal trajectory reached the connectivity threshold.');
outputT=(0:P.outputDt:P.T).';
t=unique([(0:P.denseDt:P.T).';sol.x(:);outputT]);
Z=deval(sol,t).';lambda2=zeros(size(t));
for j=1:numel(t),lambda2(j)=spectrum(P,Z(j,:).');end
[minimum,index]=min(lambda2);minimumTime=t(index);
% Refine the lowest sampled interior neighborhood on the ODE interpolant.
if index>1&&index<numel(t)
 [tm,vm]=fminbnd(@(tt)spectrum(P,deval(sol,tt)),t(index-1),t(index+1),...
  optimset('Display','off','TolX',1e-10));
 if vm<minimum,minimum=vm;minimumTime=tm;end
end
assert(minimum>P.lambdaLower);
[present,indices]=ismember(outputT,t);assert(all(present));
S=struct('t',outputT,'z',Z(indices,:),'lambda2',lambda2(indices));
R=struct('solver','ode15s','horizon',P.T,'relativeTolerance',P.RelTol,...
 'absoluteTolerance',P.AbsTol,'maximumStep',P.maxStep,'outputInterval',P.outputDt,...
 'denseInterval',P.denseDt,'acceptedNodes',numel(sol.x),'diagnosticSamples',numel(t),...
 'minimumLambda2',minimum,'minimumTime',minimumTime,'initialLambda2',spectrum(P,z0),...
 'terminalLambda2',spectrum(P,sol.y(:,end)),'terminalErrorNorm',norm(sol.y(:,end)-zs),...
 'connectivityThreshold',P.lambdaLower,'thresholdCrossings',numel(sol.xe));
writetable(table(t,lambda2,'VariableNames',{'time','lambda2'}),fullfile(folder,'lambda2.csv'));
% Strip callback handles before saving the reusable numerical interpolant.
sol.extdata=struct();save(fullfile(folder,'trajectory.mat'),'P','S','sol','-v7');
write_json(fullfile(folder,'summary.json'),R);
fprintf('%s: minimum lambda2=%.15g at t=%.9g; %d diagnostics\n',name,minimum,minimumTime,numel(t));
end

function value=spectrum(P,z)
X=reshape(z(1:50),10,5).';Y=zeros(5,2);
for i=1:5,Y(i,:)=X(i,2*i-1:2*i);end
[~,L]=audit_core('graph',P,Y);e=sort(eig((L+L.')/2));value=e(2);
end

function [value,terminal,direction]=lower_event(P,z)
value=spectrum(P,z)-P.lambdaLower;terminal=1;direction=-1;
end
