function result = reproduce_paper(outDir)
%REPRODUCE_PAPER Standalone reproduction of the manuscript simulations.
%   REPRODUCE_PAPER runs the tightened and untightened flows from scratch
%   and writes seven figures and numerical data to results/ beside this file.
%   RESULT = REPRODUCE_PAPER(OUTDIR) chooses a different output directory.
%
%   Tested environment: MATLAB R2024a with Optimization Toolbox.
%   No other user-supplied .m, .mat, data, or figure files are required.
%   Model parameters and numerical equations are preserved from:
%     simulation_with_c.m, regional_bounds.m, plot_with_c.m,
%     write_numerical_audit.m, verify_editable_figures.m.
%   Physical time is 0:0.1:1200 s. The local error tail is selected from
%   the same tightened trajectory after the original local-entry checks.
%   quadprog initializes the filter and performs offline diagnostics;
%   the integrated online controller uses the primal-dual differential flow.
%   Sampled safety checks are recorded separately from uniform certificates.
%
%   Usage: place this file in an empty folder, open MATLAB there, and type:
%       result = reproduce_paper;
%   Or open this file in the MATLAB Editor and click Run.

if nargin < 1 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
assert(exist('quadprog','file') == 2 && license('test','Optimization_Toolbox'), ...
    'MATLAB with Optimization Toolbox (quadprog) is required.');
if ~exist(outDir,'dir'), mkdir(outDir); end
runTimer = tic;
fprintf('Standalone paper reproduction: %s\n', version);
fprintf('Output directory: %s\n', outDir);
implementation = simulation_with_c(outDir, 'verify');
initial = simulation_with_c(outDir, 'initial');
assert(min([initial.hC, initial.psiC, initial.minPair, initial.minObs]) > 0, ...
    'The specified initial state failed the initial safety checks.');
result = simulation_with_c(outDir, 'full');
result.initialExtendedSafety = initial;
save(fullfile(outDir,'experiment_data.mat'), 'result', '-v7');
figureChecks = verify_editable_figures(outDir);
assert(numel(figureChecks) == 7 && all([figureChecks.editSaveReopenPass]), ...
    'Expected seven editable figure files.');
figureSummary = jsondecode(fileread(fullfile(outDir,'figure_data_audit.json')));
assert(figureSummary.entry.tailStart == 850 && figureSummary.entry.tailEnd == 1050, ...
    'The reproduced local-error window differs from the manuscript [850,1050].');
runtime = struct('matlabVersion',version,'release',version('-release'), ...
    'platform',computer,'optimizationToolbox',ver('optim'), ...
    'implementation',implementation,'initialExtendedSafety',initial, ...
    'figureCount',numel(figureChecks),'wallSeconds',toc(runTimer), ...
    'fromScratch',true,'externalUserFilesRequired',false, ...
    'status','REPRODUCTION_PASS');
fid = fopen(fullfile(outDir,'runtime_validation.json'),'w');
assert(fid >= 0, 'Cannot write runtime_validation.json.');
fileCleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(runtime,'PrettyPrint',true));
fprintf('REPRODUCTION_PASS: seven figures; tail=[850,1050]; elapsed=%.1f s\n',toc(runTimer));
end

% ======================================================================
% Embedded source: simulation_with_c.m
% ======================================================================
function result = simulation_with_c(outDir, mode, overrides)
% One regional initial condition, one controller, one continuous run per variant.
% All proposed-method plots, including the local tail, use Experiment A.
% The untightened ablation uses the same initial state and zero tightening.
% Requires MATLAB R2024a and Optimization Toolbox.
if nargin<1,outDir=fullfile(fileparts(mfilename('fullpath')),'results');end
if nargin<2,mode='full';end
if nargin<3,overrides=struct();end
if ~exist(outDir,'dir'),mkdir(outDir);end
P=parameters();fields=fieldnames(overrides);for j=1:numel(fields),P.(fields{j})=overrides.(fields{j});end
P.alpha1=(P.K2-sqrt(P.K2^2-4*P.K1))/2;P.alpha2=P.K2-P.alpha1;
lt=sqrt(1.12^2+4*.03^2);ky=P.c*P.lambdaProof/2-lt-(1.15+lt)^2/2;
assert(ky>0,'Theorem 1 consensus-gain condition failed.');
P.k=sqrt(P.c*(P.N-1)+1.15^2+lt^2/(2*ky)); % proof parameter from the NEW manuscript
C=certificates(P);dataFile=fullfile(outDir,'experiment_data.mat');
if strcmp(mode,'screen')
    result=struct('P',P,'C',C);result.untightened=simulate(P,false);
    result.metrics=screenMetrics(P,result.untightened);disp(result.metrics);
    save(dataFile,'result','-v7');return;
end
if strcmp(mode,'initial')
    X=repmat(reshape(P.y0.',1,10),5,1);[A,b,un,a]=rows(P,X,P.v0,true);
    result=struct('hC',a.hC,'psiC',sum(a.gc.*P.v0,'all')+P.alpha1*a.hC,'lambda2',a.lambda2,'minPair',min(a.hInt),'minObs',min(a.hObs));return;
end
if strcmp(mode,'certificate'),result=struct('P',P,'C',C);disp(C);save(dataFile,'result','-v7');return;end
if strcmp(mode,'verify'),result=verifyImplementation(P);disp(result);return;end
if strcmp(mode,'figures'),tmp=load(dataFile,'result');result=tmp.result;makeFigures(result,outDir);return;end
assert(any(strcmp(mode,{'full','regional','regional_refine','probe'})),'Unsupported run mode.');
result=struct('P',P,'C',C);
result.C=C;
fprintf('Shared controller: c=%g, gamma=%g, epsilon=%g, xi=%g, K1=%g, K2=%g, DeltaC=%g\n',P.c,P.gamma,P.epsilon,P.xi,P.K1,P.K2,P.DeltaC);
disp(C);save(fullfile(outDir,'parameters.mat'),'P','C');
if any(strcmp(mode,{'full','regional','regional_refine','probe'}))
    PA=P;
    if strcmp(mode,'probe'),PA.T=450;end
    if strcmp(mode,'regional_refine'),PA.RelTol=P.RelTol/10;PA.AbsTol=P.AbsTol/10;PA.maxStep=P.maxStep/2;end
    result.integrationSettings=struct('RelTol',PA.RelTol,'AbsTol',PA.AbsTol,'maxStep',PA.maxStep,'outputDt',PA.outputDt,'auditDt',PA.auditDt);
    result.proposed=simulate(PA,true);result.auditProposed=audit(PA,result.proposed,true,C);
    save(dataFile,'result','-v7');
    if true
        result.untightened=simulate(PA,false);result.auditUntightened=audit(PA,result.untightened,false,C);
        save(dataFile,'result','-v7');
    end
end
writeAudit(result,outDir);
if strcmp(mode,'full') && isfield(result,'proposed') && isfield(result,'untightened')
    assert(result.auditProposed.allQPFeasible && result.auditUntightened.allQPFeasible,'A sampled QP failed.');
    assert(result.auditProposed.sampledTighteningCheck && result.auditProposed.regionContainmentSampled,'The tightened trajectory diagnostics failed.');
    assert(result.auditProposed.designMarginInequalities && result.auditProposed.sampledBudgetCheck,'Regional design budgets failed their sampled check.');
    makeFigures(result,outDir);
end
end

function report=verifyImplementation(P)
rng(19);maxGrad=0;maxHess=0;maxPoly=0;
for k=1:8
    Y=P.ystar+.05*randn(5,2);V=randn(5,2);V=V/norm(V,'fro');[h,g,hd]=treeBarrier(P,Y,V);
    t=1e-4;hp=treeBarrier(P,Y+t*V,V);hm=treeBarrier(P,Y-t*V,V);
    maxGrad=max(maxGrad,abs((hp-hm)/(2*t)-sum(g.*V,'all')));
    maxHess=max(maxHess,abs((hp-2*h+hm)/t^2-hd));
    M=randn(4);D=randn(4);D2=randn(4);[a,b,c,cof]=detDerivatives(M,D,D2); %#ok<ASGLU>
    maxPoly=max(maxPoly,abs(a-det(M)));
end
% Independent equation-level check of c in both physical and estimator drift.
X=repmat(reshape(P.ystar.',1,10),5,1)+.02*randn(5,10);V=.03*randn(5,2);
[~,~,~,flow]=rows(P,X,V,true);[~,L]=graph(P,flow.Y);
expectedX=zeros(5,10);expectedDrift=-P.gamma*V;
for i=1:5
    consensus=zeros(1,10);
    for j=1:5,consensus=consensus+P.c*flow.W(i,j)*(X(i,:)-X(j,:));end
    expectedX(i,:)=-consensus;ii=2*i-1:2*i;expectedX(i,ii)=V(i,:);
    expectedDrift(i,:)=expectedDrift(i,:)-consensus(ii);
end
flowError=max(abs([flow.Xdot(:)-expectedX(:);flow.f0(:)-expectedDrift(:)]));
report=struct('maxGradientDirectionalError',maxGrad,'maxHessianDirectionalError',maxHess,'maxPolynomialDetError',maxPoly,...
    'consensusGainEquationError',flowError,'pass',maxGrad<1e-6 && maxHess<1e-5 && maxPoly<1e-10 && flowError<1e-11);
assert(report.pass,'Derivative verification failed.');
end

function P=parameters()
P.N=5;P.n=2;P.edges=nchoosek(1:5,2);P.m=10;
P.y0=[-3 .002;-3 -1.15;-3 -.55;-3 .55;-3 1.15];P.v0=zeros(5,2);
P.ystar=[3 0;3 -.55;3 -1.15;3 1.15;3 .55];
P.obs=[.24 0;1.72 1.1;1.72 -1.1];P.rObs=[.625;.315;.315];
P.dsafe=.345;P.dc=3.2;P.lambdaLower=.2;
P.alpha=.03;P.B=1.15*eye(5)-.03*ones(5);P.p=P.B*P.ystar;
P.c=40;P.gamma=45;P.lambdaProof=.2;P.lambdaLocal=3;P.xi=1.4;P.epsilon=2e-5;
P.K1=14.5;P.K2=8.8;P.alpha1=(P.K2-sqrt(P.K2^2-4*P.K1))/2;P.alpha2=P.K2-P.alpha1;
P.DeltaInt=1.2;P.DeltaObs=1.8;P.DeltaC=2.5;P.barrierScale=.06;P.connScale=.1; % physical margin equivalents
P.nu=.02;P.theta=50;P.deltaP=.001;P.localRadius=2e-4;
P.T=1200;P.outputDt=.1;P.maxStep=.25;P.tailStart=800;P.tailDuration=200;
P.RelTol=2e-8;P.AbsTol=1e-11;
P.auditDt=.1;P.regionDiameter=2.5;P.regionBox=[-3.1 3.1;-1.3 1.3];P.regionSpeedMax=.5;P.auditInputRadius=15;
P.trackingBudget=.165;P.connMismatchBudget=.055; % design budgets, not proven uniform bounds
P.nU=10;P.nZ=25;P.nLambda=45;P.dim=140;
P.iX=1:50;P.iV=51:60;P.iU=61:70;P.iZ=71:95;P.iL=96:140;
P.treeScale=512;
end

function C=certificates(P)
N=P.N;C.mu=1;C.l=1.15;C.tildeL=sqrt(1.12^2+4*.03^2);
C.cLower=(2*C.tildeL+(C.l+C.tildeL)^2/C.mu)/P.lambdaProof;
C.kappaY=P.c*P.lambdaProof/2-C.tildeL-(C.l+C.tildeL)^2/(2*C.mu);
C.kProof=sqrt(P.c*(N-1)+C.l^2/C.mu+C.tildeL^2/(2*C.kappaY));
C.gammaLower=2*C.kProof;
C.nominalCoefficients=[P.k*C.mu/(4*N),P.k*C.kappaY/2,P.gamma-2*P.k];
C.M=diag(C.nominalCoefficients);C.alphaNom=min(C.nominalCoefficients);
Pk=[P.k*P.gamma P.k;P.k 1];
C.mUnder=.5*min([eig(Pk);P.k]);C.mBar=.5*max([eig(Pk);P.k]);
C.nominalDisplayedGainConditions=P.c>C.cLower && P.gamma>C.gammaLower && C.alphaNom>0;
% Independent quadratic-game cross-check of the new general nominal bound.
R=zeros(5,25);D=R;for i=1:5,R(i,(i-1)*5+i)=1;D(i,(i-1)*5+(1:5))=-.03;D(i,(i-1)*5+i)=1.12;end
Qx=P.k*P.c/2*P.lambdaLocal*kron(eye(5)-ones(5)/5,eye(5))+P.k/2*(R'*D+D'*R)-D'*D;
C.exactGameAlpha=min([eig(Qx);P.gamma-P.k-P.c*5/(2*P.k)-.25]);
C.Lz=sqrt((1+8*P.c)^2+(P.gamma+8*P.c+C.tildeL)^2);C.Lu=1;C.cf=.5;
C.sigma=(1+P.k)/(C.tildeL*C.Lz);
C.epsUpper=4*C.cf*C.alphaNom*C.sigma/((1+P.k+C.sigma*C.tildeL*C.Lz)^2+4*C.alphaNom*C.sigma*C.tildeL);
C.Q=[C.alphaNom,-.5*(1+P.k+C.sigma*C.tildeL*C.Lz);...
 -.5*(1+P.k+C.sigma*C.tildeL*C.Lz),C.sigma*(C.cf/P.epsilon-C.tildeL)];
C.minEigQ=min(eig(C.Q));
assert(C.nominalDisplayedGainConditions && P.epsilon<C.epsUpper && C.minEigQ>0,'New theorem parameter conditions failed.');
X=repmat(reshape(P.ystar.',1,10),5,1);V=zeros(5,2);
[A,b,nom,aux]=rows(P,X,V,true);C.targetMinReducedSlack=min(-b);C.targetU=norm(nom);
C.targetTreeMargin=aux.hC;C.targetLambda2=aux.lambda2;
% Global derivative bounds for h=tree/512-lambdaLower, using 125 trees.
Gedge=8/(exp(1)*P.dc);
Hedge=(8/exp(1)+4*(256/exp(3)+54/exp(2)))/P.dc^2;
C.GhLoose=625/512*4*sqrt(2)*Gedge;
% For each spanning tree, ||B_tree||<=sqrt(N). With s_e=d_e^2/(dc^2-d_e^2),
% sum s_e*(1+s_e)^3 <= S*(1+S)^3, S=sum s_e. Max occurs at S=(1+sqrt(3))/2.
sm=(1+sqrt(3))/2;C.Gh=625/512*(2*sqrt(5)/P.dc)*sqrt(sm*(1+sm)^3)*exp(-sm);
C.Hh=625/512*(8*Hedge+24*Gedge^2);
C.rTheta=P.localRadius;
C.localGraphPerturbationBound=2*P.m*sqrt(2)*Gedge*C.rTheta;
C.localGraphLower=C.targetLambda2-C.localGraphPerturbationBound;
assert(C.localGraphLower>P.lambdaLocal,'Local graph lower bound failed.');
r=C.rTheta;dmax=max(vecnorm(P.ystar(P.edges(:,1),:)-P.ystar(P.edges(:,2),:),2,2))+sqrt(2)*r;
hmin=min(aux.hInt);obsDist=max(vecnorm(reshape(P.ystar,[5,1,2])-reshape(P.obs,[1,3,2]),2,3),[],'all')+r;
% Conservative unweighted local pair/obstacle slack bounds; both split rows.
pairSlack=.5*(P.K1*(hmin-2*dmax*sqrt(2)*r)-P.DeltaInt)...
 -2*dmax*(P.gamma+8*P.c+C.tildeL+P.K2)*r-2*r^2;
chiMin=exp(-dmax^2/(P.dc^2-dmax^2));
obsSlack=P.K1*(min(aux.hObs)-2*obsDist*r)-P.DeltaObs...
 -2*obsDist*(P.gamma+8*P.c+C.tildeL+P.K2)*r-2*r^2;
connSlack=(P.K1*(aux.hC-C.Gh*r)-P.DeltaC)/5 ...
 -C.Gh*(P.gamma+8*P.c+C.tildeL+P.K2)*r-C.Hh*((1+8*P.c)*r)^2/5;
C.etaI=min([P.barrierScale*[chiMin*pairSlack,obsSlack],P.connScale*connSlack,P.deltaP]);
C.localBarrierLower=[hmin-2*dmax*sqrt(2)*r,min(aux.hObs)-2*obsDist*r,aux.hC-C.Gh*r];
C.localExtendedLower=P.alpha1*C.localBarrierLower-[2*dmax*sqrt(2)*r,2*obsDist*r,C.Gh*r];
assert(min(C.localExtendedLower)>0,'Local ball not contained in extended safe set.');
% Row-norm bounds over the local ball, including mismatch coefficients.
pairRow=sqrt((P.barrierScale*2*dmax)^2+2);obsRow=P.barrierScale*2*obsDist;
connRow=sqrt((P.connScale*C.Gh)^2+8^2);passRow=0;
C.Bbar=sqrt(20*pairRow^2+15*obsRow^2+5*connRow^2+5*passRow^2);
C.rFast=1e-3;
% Rowwise bounds retain the zero coefficient of the passive rows.
% This proves exactly g<=-etaI/2 required by the manuscript's local tube.
C.localNominalSlackLower=[P.barrierScale*[chiMin*pairSlack,obsSlack],P.connScale*connSlack,P.deltaP];
C.localRowNormUpper=[pairRow,obsRow,connRow,passRow];
C.rFast=min(1e-3,.5*min((C.localNominalSlackLower(1:3)-C.etaI/2)./C.localRowNormUpper(1:3)));
C.fastTubeSlackLower=C.localNominalSlackLower-C.localRowNormUpper*C.rFast;
assert(min(C.fastTubeSlackLower)>=C.etaI/2,'Local strict inactivity tube failed.');
C.kappa=2+(C.Bbar^2+1)*C.rFast/C.etaI;
C.cfRecomputed=min(.5,C.kappa*C.etaI/(2*C.rFast)-C.Bbar^2/2);
C.mUnderComposite=min(C.mUnder,.5*C.sigma*min(1,C.kappa));
C.mBarComposite=max(C.mBar,.5*C.sigma*max(1,C.kappa));
C.rhoComposite=.25*min(C.mUnder*r^2,.5*C.sigma*min(1,C.kappa)*C.rFast^2);
C.rateBound=C.minEigQ/(2*C.mBarComposite);
C.localAnalyticConditions=C.nominalDisplayedGainConditions && C.exactGameAlpha>=C.alphaNom && C.etaI>0 && C.minEigQ>0 && r<=P.nu/2 && C.cfRecomputed>=C.cf;
assert(C.localAnalyticConditions,'Local theorem certificate failed.');
C.regionalUniformCertificate=false;
C.notes='Nominal gains and Lyapunov coefficients use the revised manuscript with consensus gain c and lambdaProof=0.2. Local strict-inactivity bounds use lambdaLocal=3. Regional uniform feasibility/tracking are not established by sampled maxima. Corollary 1 is not invoked.';
end

function res=simulate(P,tight)
% Integrate each smooth face and locate changes of the projected dual flow
% with events. This avoids finite-difference chattering at lambda=0.
X0=repmat(reshape(P.y0.',1,10),5,1);V0=P.v0;
[A,b,nom]=rows(P,X0,V0,tight);H=diag([ones(10,1);P.xi*ones(25,1)]);
opt=optimoptions('quadprog','Display','off','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
[w,exitflag,lm]=solveScaledQP(P,A,b,nom);
assert(exitflag>0,'Initial QP infeasible.');la=max(0,lm.ineqlin);la(A*w+b<-1e-7)=0;
zz0=[reshape(X0.',50,1);reshape(V0.',10,1);w;la];active=la>1e-9;zz0(P.iL(~active))=0;
gridT=(0:P.outputDt:P.T).';allZ=nan(numel(gridT),140);allZ(1,:)=zz0.';
t0=0;start=tic;lastPrint=-100;lastRhs=0;nEvents=0;eventLog=zeros(0,2);
fprintf('Integrating projected flow by smooth faces: tight=%d\n',tight);
fprintf('Initial filter distortion %.9g, active multipliers %d\n',norm(w(1:10)-nom),nnz(active));
while t0<P.T
    opts=odeset('RelTol',P.RelTol,'AbsTol',P.AbsTol,'MaxStep',P.maxStep,...
        'Jacobian',@jacobian,'Events',@events,'OutputFcn',@progress);
    sol=ode15s(@rhs,[t0 P.T],zz0,opts);tend=sol.x(end);
    pick=find(gridT>=t0 & gridT<=tend);if ~isempty(pick),allZ(pick,:)=deval(sol,gridT(pick)).';end
    if tend>=P.T,break;end
    if isempty(sol.ie),error('Solver terminated without reaching horizon or a face event.');end
    hit=unique(sol.ie(sol.xe==sol.xe(end)));zz0=sol.y(:,end);zz0(P.iL(hit))=0;
    active(hit)=~active(hit);zz0(P.iL(~active))=0;
    nEvents=nEvents+numel(hit);eventLog=[eventLog;[repmat(tend,numel(hit),1),hit(:)]]; %#ok<AGROW>
    if mod(nEvents,20)==0,fprintf('  event %d at t=%.9g, active=%d\n',nEvents,tend,nnz(active));end
    assert(nEvents<3000,'Excessive face switching; inspect integration tolerances.');
    t0=tend;
end
assert(all(isfinite(allZ),'all'),'Missing trajectory output samples.');
res=struct('t',gridT,'z',allZ,'tight',tight,'wallSeconds',toc(start),'initialQPExit',exitflag,...
    'completed',true,'faceEvents',eventLog,'minDualNumerical',min(allZ(:,P.iL),[],'all'));
res.controller=controller(P);if ~tight,res.controller.DeltaInt=0;res.controller.DeltaObs=0;res.controller.DeltaC=0;end
fprintf('Finished tight=%d, t=%.0f, events=%d, wall=%.1f s\n',tight,P.T,nEvents,toc(start));
    function dz=rhs(tt,zz)
        if toc(start)-lastRhs>30
            fprintf('  internal t=%.9g, active=%d, maxLambda=%.4g\n',tt,nnz(active),max(zz(P.iL)));lastRhs=toc(start);
        end
        X=reshape(zz(P.iX),10,5).';V=reshape(zz(P.iV),2,5).';[AA,bb,un,aux]=rows(P,X,V,tight);
        ww=zz([P.iU P.iZ]);laa=zz(P.iL);laa(~active)=0;
        dw=-H*ww+[un;zeros(25,1)]-AA.'*laa;dl=(AA*ww+bb).*active;
        acc=aux.f0+reshape(ww(1:10),2,5).';
        dz=[reshape(aux.Xdot.',50,1);reshape(acc.',10,1);dw/P.epsilon;dl/P.epsilon];
    end
    function J=jacobian(tt,zz)
        J=zeros(140);X=reshape(zz(P.iX),10,5).';V=reshape(zz(P.iV),2,5).';[AA,~]=rows(P,X,V,tight);
        for jj=1:60
            hh=1e-6*(1+abs(zz(jj)));zp=zz;zm=zz;zp(jj)=zp(jj)+hh;zm(jj)=zm(jj)-hh;
            J(:,jj)=(rhs(tt,zp)-rhs(tt,zm))/(2*hh);
        end
        J(P.iV,P.iU)=eye(10);J([P.iU P.iZ],[P.iU P.iZ])=-H/P.epsilon;
        J([P.iU P.iZ],P.iL)=-AA.'/P.epsilon;J([P.iU P.iZ],P.iL(~active))=0;
        J(P.iL,[P.iU P.iZ])=(AA.*active)/P.epsilon;
    end
    function [value,isterminal,direction]=events(tt,zz)
        X=reshape(zz(P.iX),10,5).';V=reshape(zz(P.iV),2,5).';[AA,bb]=rows(P,X,V,tight);
        value=AA*zz([P.iU P.iZ])+bb;value(active)=zz(P.iL(active));
        direction=ones(45,1);direction(active)=-1;isterminal=ones(45,1);
        if tt<=t0+4*eps(max(1,t0))
            value(active)=max(value(active),1e-12);value(~active)=min(value(~active),-1e-12);
        end
    end
    function stop=progress(tt,~,flag)
        stop=false;if isempty(flag)&&~isempty(tt)&&tt(end)-lastPrint>=100
            fprintf('  tight=%d t=%.3f / %.0f, wall=%.1f s\n',tight,tt(end),P.T,toc(start));lastPrint=tt(end);
        end
    end
end

function [A,b,nom,aux]=rows(P,X,V,tight)
Y=zeros(5,2);for i=1:5,Y(i,:)=X(i,2*i-1:2*i);end
[W,L]=graph(P,Y);rho=P.c*L*X;Xdot=-rho;
f0=-P.gamma*V;grad=zeros(5,2);
for i=1:5
    ii=2*i-1:2*i;Xdot(i,ii)=V(i,:);f0(i,:)=f0(i,:)-rho(i,ii);
    localY=reshape(X(i,:),2,5).';grad(i,:)=P.B(i,:)*localY-P.p(i,:);
end
nom=-reshape(grad.',10,1);A=zeros(45,35);b=zeros(45,1);hInt=zeros(10,1);hObs=zeros(15,1);rr=0;
for e=1:10
    i=P.edges(e,1);j=P.edges(e,2);d=Y(i,:)-Y(j,:);v=V(i,:)-V(j,:);h=sum(d.^2)-P.dsafe^2;hInt(e)=h;
    for side=1:2
        if side==1,a=i;sgn=1;other=2*e;own=2*e-1;else,a=j;sgn=-1;other=2*e-1;own=2*e;end
        g=sgn*2*d;rr=rr+1;A(rr,2*a-1:2*a)=-W(i,j)*g;
        A(rr,10+own)=W(i,j);A(rr,10+other)=-W(i,j);
        b(rr)=W(i,j)*(-g*f0(a,:).'-P.K2*g*V(a,:).'-sum(v.^2)-.5*P.K1*h+.5*tight*P.DeltaInt);
    end
end
for i=1:5
    for o=1:3
        d=Y(i,:)-P.obs(o,:);h=sum(d.^2)-P.rObs(o)^2;hObs((i-1)*3+o)=h;g=2*d;
        rr=rr+1;A(rr,2*i-1:2*i)=-g;
        b(rr)=-g*f0(i,:).'-2*sum(V(i,:).^2)-P.K2*g*V(i,:).'-P.K1*h+tight*P.DeltaObs;
    end
end
for i=1:5
    Yi=reshape(X(i,:),2,5).';Vi=reshape(Xdot(i,:),2,5).';
    [hc,gc,Hdir]=treeBarrier(P,Yi,Vi);
    rr=rr+1;A(rr,2*i-1:2*i)=-gc(i,:);A(rr,31:35)=L(i,:);
    b(rr)=-gc(i,:)*f0(i,:).'-P.K2*gc(i,:)*V(i,:).'-(Hdir+P.K1*hc)/5+tight*P.DeltaC/5;
end
% Positive barrier normalization: pair/obstacle use barrierScale, connectivity
% uses connScale. Each Delta has the same scaling as its barrier.
% Mismatch matrices retain the manuscript form in both experiments/ablations.
A(1:35,1:10)=P.barrierScale*A(1:35,1:10);b(1:35)=P.barrierScale*b(1:35);
A(36:40,1:10)=P.connScale*A(36:40,1:10);b(36:40)=P.connScale*b(36:40);
for i=1:5
    rr=rr+1;s=sum(V(i,:).^2);tau=max(0,min(1,(4*s-P.nu^2)/(3*P.nu^2)));
    chi=10*tau^3-15*tau^4+6*tau^5;ii=2*i-1:2*i;
    A(rr,ii)=chi*V(i,:);b(rr)=-chi*(V(i,:)*nom(ii)+P.theta*s)-(1-chi)*P.deltaP;
end
if nargout>=4
    [hc,gc,Hdir]=treeBarrier(P,Y,V);ev=sort(eig(L));
    aux=struct('Y',Y,'f0',f0,'Xdot',Xdot,'hInt',hInt,'hObs',hObs,'hC',hc,'lambda2',ev(2),...
        'gc',gc,'Hdir',Hdir,'W',W);
end
end

function [W,L]=graph(P,Y)
D=P.dc^2;dx=Y(:,1)-Y(:,1).';dy=Y(:,2)-Y(:,2).';s=dx.^2+dy.^2;
W=zeros(5);m=s<D;W(m)=exp(-s(m)./(D-s(m)));W(1:6:end)=0;L=diag(sum(W,2))-W;
end

function [h,g,Hdir]=treeBarrier(P,Y,V)
% Exact polynomial determinant derivatives: no eigenvector differentiation,
% finite-difference Hessians or derivative clipping.
D=P.dc^2;E=P.edges;inc=zeros(5,10);w=zeros(10,1);wd=w;wdd=w;gs=zeros(10,2);
for e=1:10
    i=E(e,1);j=E(e,2);inc(i,e)=1;inc(j,e)=-1;d=Y(i,:)-Y(j,:);v=V(i,:)-V(j,:);s=sum(d.^2);
    if s<D
        w(e)=exp(-s/(D-s));ws=-D/(D-s)^2*w(e);wss=(D^2/(D-s)^4-2*D/(D-s)^3)*w(e);
        sd=2*d*v.';wd(e)=ws*sd;wdd(e)=wss*sd^2+2*ws*sum(v.^2);gs(e,:)=2*ws*d;
    end
end
J=inc(1:4,:);M=(J.*w.')*J.';Md=(J.*wd.')*J.';Mdd=(J.*wdd.')*J.';
detM=det(M);h=5*detM/P.treeScale-P.lambdaLower;g=zeros(5,2);
if rcond(M)>1e-11
    Z=M\Md;val=5*detM/P.treeScale;
    Hdir=val*(trace(Z)^2+trace(M\Mdd)-trace(Z*Z));
    R=M\J;coeff=sum(J.*R,1)*val;
else
    [val,d1,d2,cof]=detDerivatives(M,Md,Mdd); %#ok<ASGLU>
    Hdir=5*d2/P.treeScale;coeff=zeros(1,10);
    for e=1:10,coeff(e)=5/P.treeScale*sum(cof.*(J(:,e)*J(:,e).'),'all');end
end
for e=1:10,g(E(e,1),:)=g(E(e,1),:)+coeff(e)*gs(e,:);g(E(e,2),:)=g(E(e,2),:)-coeff(e)*gs(e,:);end
end

function [val,d1,d2,cof]=detDerivatives(M,D,D2)
p=perms(1:4);val=0;d1=0;d2=0;cof=zeros(4);
for k=1:24
    a=p(k,:);sgn=(-1)^sum(triu(a.'>a,1),'all');idx=sub2ind([4 4],1:4,a);m=M(idx);d=D(idx);dd=D2(idx);
    val=val+sgn*prod(m);
    for i=1:4
        rem=setdiff(1:4,i);c=sgn*prod(m(rem));cof(i,a(i))=cof(i,a(i))+c;d1=d1+c*d(i);d2=d2+c*dd(i);
        for j=rem,d2=d2+sgn*d(i)*d(j)*prod(m(setdiff(rem,j)));end
    end
end
end

function audit=audit(P,res,tight,C)
stride=max(1,round(P.auditDt/P.outputDt));ind=unique([1:stride:numel(res.t),numel(res.t)]);ns=numel(ind);
ru=zeros(ns,1);distort=ru;deltaC=ru;deltaCBall=ru;uNorm=ru;uStarNorm=ru;minG=ru;feas=ru;Vnom=ru;Vc=ru;err=ru;hI=ru;hO=ru;hC=ru;lam2=ru;resQP=ru;active=ru;maxTrueQ=ru;stationarity=ru;complementarity=ru;
H=diag([ones(10,1);P.xi*ones(25,1)]);opt=optimoptions('quadprog','Display','off','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
fprintf('Offline diagnostic QPs: tight=%d, %d samples\n',tight,ns);
for a=1:ns
    if mod(a-1,2000)==0,fprintf('  audit sample %d / %d\n',a,ns);end
    z=res.z(ind(a),:).';X=reshape(z(P.iX),10,5).';V=reshape(z(P.iV),2,5).';
    [A,b,un,aux]=rows(P,X,V,tight);u=z(P.iU);w=z([P.iU P.iZ]);la=z(P.iL);
    [ws,exitflag,lm]=solveScaledQP(P,A,b,un);feas(a)=exitflag;
    stationarity(a)=norm(H*ws-[un;zeros(25,1)]+A.'*lm.ineqlin,inf);
    complementarity(a)=norm(lm.ineqlin.*(A*ws+b),inf);
    if exitflag>0,ru(a)=norm(u-ws(1:10));distort(a)=norm(ws(1:10)-un);resQP(a)=max(A*ws+b);uStarNorm(a)=norm(ws(1:10));else,ru(a)=NaN;distort(a)=NaN;resQP(a)=NaN;end
    trueConn=-sum(aux.gc.*(aux.f0+reshape(u,2,5).'),'all')-aux.Hdir-P.K2*sum(aux.gc.*V,'all')-P.K1*aux.hC;
    hatConn=sum(A(36:40,1:10)*u+b(36:40))/P.connScale-tight*P.DeltaC;deltaC(a)=abs(trueConn-hatConn);
    da=-reshape(aux.gc.',1,10)-sum(A(36:40,1:10),1)/P.connScale;
    db=trueConn-hatConn-da*u;deltaCBall(a)=abs(db)+P.auditInputRadius*norm(da);uNorm(a)=norm(u);
    minG(a)=max(A*w+b);active(a)=sum(la>1e-6);
    eX=X-repmat(reshape(P.ystar.',1,10),5,1);eY=aux.Y-P.ystar;
    Vnom(a)=.5*sum(V.^2,'all')+P.k*sum(eY.*V,'all')+P.k*P.gamma/2*sum(eY.^2,'all')+P.k/2*(sum(eX.^2,'all')-sum(eY.^2,'all'));
    Vf=.5*norm([u-un;z(P.iZ)])^2+C.kappa/2*norm(la)^2;Vc(a)=Vnom(a)+C.sigma*Vf;
    err(a)=norm([eX(:);V(:);u-un;z(P.iZ);la]);
    hI(a)=min(aux.hInt);hO(a)=min(aux.hObs);hC(a)=aux.hC;lam2(a)=aux.lambda2;
    qobs=(A(21:35,1:10)*u+b(21:35))/P.barrierScale-tight*P.DeltaObs;
    qpair=zeros(10,1);for e=1:10,rr=2*e-1:2*e;i=P.edges(e,1);j=P.edges(e,2);qpair(e)=(sum(A(rr,1:10)*u+b(rr))/(P.barrierScale*max(aux.W(i,j),realmin)))-tight*P.DeltaInt;end
    maxTrueQ(a)=max([qobs;qpair;trueConn]);
end
audit=struct('samplePeriod',P.outputDt*stride,'allQPFeasible',all(feas>0),'maxTrackingSampled',max(ru),...
 'maxDistortionSampled',max(distort),'deltaCSampled',max(deltaC),'minInterSampled',min(hI),'minObstacleSampled',min(hO),...
 'minTreeSampled',min(hC),'minLambda2Sampled',min(lam2),'maxTrueECBFResidualSampled',max(maxTrueQ),'maxLocalRowResidualSampled',max(minG),...
 'maxQPResidual',max(resQP),'maxActiveMultipliers',max(active),'finalError',err(end),'regionalUniformCertificate',false);
audit.maxQPStationarityResidual=max(stationarity);audit.maxQPComplementarityResidual=max(complementarity);
audit.deltaCForInputBallOverSampledStates=max(deltaCBall);audit.inputBallRadius=P.auditInputRadius;audit.maxInputNormSampled=max(uNorm);audit.maxOptimizerInputNormSampled=max(uStarNorm);
audit=regional_bounds(P,res,C,audit);
audit.t=res.t(ind);audit.ru=ru;audit.distortion=distort;audit.deltaC=deltaC;audit.Vnom=Vnom;audit.Vc=Vc;audit.error=err;
audit.hInt=hI;audit.hObs=hO;audit.hTree=hC;audit.lambda2=lam2;
audit.corollary1Claimed=false;
audit.certificateScope='Pointwise sampled diagnostics; these maxima are not uniform regional certificates.';
disp(rmfield(audit,{'t','ru','distortion','deltaC','Vnom','Vc','error','hInt','hObs','hTree','lambda2'}));
end

function makeFigures(R,outDir)
plot_with_c(R,outDir);
end

function S=controller(P)
fields={'c','gamma','epsilon','xi','K1','K2','DeltaInt','DeltaObs','DeltaC','nu','theta','deltaP','dc','dsafe','lambdaLower','barrierScale','connScale'};
S=struct();for j=1:numel(fields),S.(fields{j})=P.(fields{j});end
end

function writeAudit(R,outDir)
write_numerical_audit(R,outDir);
end

function [w,flag,lm]=solveScaledQP(P,A,b,nom)
% Exact variable/constraint change for the OFFLINE QP only.
% It leaves the original objective, feasible set, solution and KKT multiplier
% mapping unchanged; the online primal-dual differential equation is untouched.
beta=sqrt(P.barrierScale);bc=sqrt(P.connScale);
T=[ones(10,1);beta*ones(20,1);bc*ones(5,1)];row=[ones(35,1)/beta;ones(5,1)/bc;ones(5,1)];
At=(A.*T.').*row;bt=b.*row;Ht=diag([ones(10,1);P.xi*beta^2*ones(20,1);P.xi*bc^2*ones(5,1)]);
opt=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex',...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000);
[x,~,flag,~,ll]=quadprog(Ht,[-nom;zeros(25,1)],At,-bt,[],[],[],[],[nom;zeros(25,1)],opt);
w=T.*x;lm=struct('ineqlin',row.*ll.ineqlin);
end

function S=screenMetrics(P,R)
nt=numel(R.t);Y=reshape(R.z(:,[1 2 13 14 25 26 37 38 49 50]),nt,2,5);Y=permute(Y,[1 3 2]);
hp=zeros(nt,10);ho=zeros(nt,5,3);hc=zeros(nt,1);la=hc;
for j=1:10,ij=P.edges(j,:);hp(:,j)=sum((Y(:,ij(1),:)-Y(:,ij(2),:)).^2,3)-P.dsafe^2;end
for i=1:5,for j=1:3,ho(:,i,j)=sum((Y(:,i,:)-reshape(P.obs(j,:),1,1,2)).^2,3)-P.rObs(j)^2;end,end
for k=1:nt,y=squeeze(Y(k,:,:));[~,L]=graph(P,y);ev=sort(eig(L));la(k)=ev(2);hc(k)=5*det(L(1:4,1:4))/P.treeScale-P.lambdaLower;end
[mc,ic]=min(hc);S=struct('minTree',mc,'treeMinTime',R.t(ic),'minLambda2',min(la),'minPair',min(hp,[],'all'),'minObs',min(ho,[],'all'),'initialTree',hc(1),'initialLambda2',la(1));
end

% ======================================================================
% Embedded source: regional_bounds.m
% ======================================================================
function A=regional_bounds(P,S,C,A)
% Analytic coefficient bounds on an explicitly declared regional geometry.
% Feasibility and tracking values in A remain sampled, not uniform proofs.
A.regionDiameter=P.regionDiameter;A.regionBox=P.regionBox;
A.barAIntBox=2*P.regionDiameter;
farthest=zeros(size(P.obs,1),1);
for j=1:size(P.obs,1)
    dx=max(abs(P.regionBox(1,:)-P.obs(j,1)));
    dy=max(abs(P.regionBox(2,:)-P.obs(j,2)));
    farthest(j)=hypot(dx,dy);
end
A.barAObsBox=2*max(farthest);A.barAC=C.Gh;
A.lhsIntWithSampledRu=sqrt(2)*A.barAIntBox*A.maxTrackingSampled;
A.lhsObsWithSampledRu=A.barAObsBox*A.maxTrackingSampled;
A.lhsConnWithSampledBounds=sqrt(P.N)*C.Gh*A.maxTrackingSampled+A.deltaCForInputBallOverSampledStates;
A.sampledTighteningCheck=S.tight && A.lhsIntWithSampledRu<=P.DeltaInt && A.lhsObsWithSampledRu<=P.DeltaObs && A.lhsConnWithSampledBounds<=P.DeltaC;
A.trackingDesignBudget=P.trackingBudget;A.connMismatchDesignBudget=P.connMismatchBudget;
A.designLHS=[sqrt(2)*A.barAIntBox*P.trackingBudget,A.barAObsBox*P.trackingBudget,sqrt(P.N)*C.Gh*P.trackingBudget+P.connMismatchBudget];
A.designRHS=[P.DeltaInt,P.DeltaObs,P.DeltaC];
A.designMarginInequalities=all(A.designLHS<A.designRHS);
A.sampledBudgetCheck=S.tight && A.maxTrackingSampled<P.trackingBudget && A.deltaCForInputBallOverSampledStates<P.connMismatchBudget;
A.safeActivationRelativeSpeedBound=2*P.regionSpeedMax;
A.safeActivationLHS=P.alpha1*(P.dc^2-P.dsafe^2);
A.safeActivationRHS=2*P.dc*A.safeActivationRelativeSpeedBound;
Y=zeros(numel(S.t),P.N,2);
for i=1:P.N,j=(i-1)*10+2*i-1;Y(:,i,:)=S.z(:,j:j+1);end
maxDiameter=0;
for e=1:size(P.edges,1)
    ij=P.edges(e,:);maxDiameter=max(maxDiameter,max(sqrt(sum((Y(:,ij(1),:)-Y(:,ij(2),:)).^2,3))));
end
V=reshape(S.z(:,51:60),[],2,5);
A.maxPairDistanceSampled=maxDiameter;A.maxAgentSpeedSampled=max(sqrt(sum(V.^2,2)),[],'all');
A.positionMinSampled=[min(Y(:,:,1),[],'all'),min(Y(:,:,2),[],'all')];
A.positionMaxSampled=[max(Y(:,:,1),[],'all'),max(Y(:,:,2),[],'all')];
A.regionContainmentSampled=maxDiameter<=P.regionDiameter && A.maxAgentSpeedSampled<=P.regionSpeedMax && ...
    all(A.positionMinSampled>=P.regionBox(:,1).') && all(A.positionMaxSampled<=P.regionBox(:,2).');
A.regionalUniformCertificate=false;
end

% ======================================================================
% Embedded source: plot_with_c.m
% ======================================================================
function report=plot_with_c(R,outDir)
% Seven figure files, exclusively from the stored regional runs.
% Figures for the revised consensus-gain-c model in this isolated package.
% With no arguments, load results/experiment_data.mat in this package.
% Trajectories: my_trajectory.fig (tightened), my_trajectory_2.fig (untightened).
% Position components: HB_flow_hebing.fig, one row and two columns.
% The error tail is selected after certified local-entry checks on R.proposed.
if nargin<1 || isempty(R)
    source=fullfile(fileparts(mfilename('fullpath')),'results','experiment_data.mat');
    assert(isfile(source),'Run RUN_ALL first, or pass R and outDir.');
    saved=load(source,'result');R=saved.result;
end
if nargin<2 || isempty(outDir),outDir=fullfile(fileparts(mfilename('fullpath')),'results');end
if ~exist(outDir,'dir'),mkdir(outDir);end
P=R.P;T=geometry(P,R.proposed);U=geometry(P,R.untightened);
keys={'c','gamma','epsilon','xi','K1','K2','nu','theta','deltaP','dc','dsafe','lambdaLower','barrierScale','connScale'};
for j=1:numel(keys),assert(isequal(R.proposed.controller.(keys{j}),R.untightened.controller.(keys{j})));end
old=get(groot,{'defaultFigureVisible','defaultAxesFontName','defaultAxesFontSize','defaultLineLineWidth'});
restore=onCleanup(@()set(groot,{'defaultFigureVisible','defaultAxesFontName','defaultAxesFontSize','defaultLineLineWidth'},old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Times New Roman','defaultAxesFontSize',10,'defaultLineLineWidth',1.4);
cols=[.85 .13 .10;0 .46 .22;.12 .32 .82;.66 .15 .63;0 .56 .64];
styles={'-','--','-.',':','-'};

% Two independent trajectory figures. Each has one obstacle-boundary inset.
% Both insets refer to agent 1 and the same central obstacle, at the closest
% saved sample of the corresponding trajectory. Axes show their actual scale.
zoomDetails=repmat(struct(),2,1);
for mode=1:2
    if mode==1
        D=T;label='Tightened';fileName='my_trajectory';
    else
        D=U;label='Untightened';fileName='my_trajectory_2';
    end
    fig=newFigure([label ' trajectories'],7.16,3.7);
    ax=axes(fig,'Position',[.08 .21 .60 .68]);hold(ax,'on');obstacles(ax,P);
    hh=gobjects(5,1);
    for i=1:5
        hh(i)=plot(ax,D.Y(:,i,1),D.Y(:,i,2),'Color',cols(i,:),'LineStyle',styles{i},'LineWidth',1.5,'DisplayName',sprintf('Agent %d',i));
        plot(ax,P.y0(i,1),P.y0(i,2),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerSize',4,'HandleVisibility','off');
        plot(ax,P.ystar(i,1),P.ystar(i,2),'d','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerSize',4,'HandleVisibility','off');
        xx=D.Y(:,i,1);yy=D.Y(:,i,2);bad=min(D.obsAll(:,i,:),[],3)<0;xx(~bad)=NaN;yy(~bad)=NaN;
        plot(ax,xx,yy,'k-','LineWidth',2.1,'HandleVisibility','off');
    end
    axis(ax,'equal');xlim(ax,[-3.25 3.25]);ylim(ax,[-1.65 1.65]);styleAxes(ax);
    xlabel(ax,'x (m)');ylabel(ax,'y (m)');title(ax,label,'FontWeight','normal');
    lg=legend(ax,hh,'Location','none','Orientation','horizontal','NumColumns',5,'Box','off','FontSize',9);
    lg.Position=[.12 .025 .78 .075];

    ia=1;io=1;xy=squeeze(D.Y(:,ia,:));oc=P.obs(io,:);
    distances=vecnorm(xy-oc,2,2);[closest,it]=min(distances);
    clearance=closest-P.rObs(io);pt=xy(it,:);
    boundary=oc+P.rObs(io)*(pt-oc)/closest;
    center=(pt+boundary)/2;half=max(.005,1.15*abs(clearance));
    if mode==1
        assert(clearance>0,'The tightened inset must have positive obstacle clearance.');
        distanceLabel=sprintf('Clearance: %.2f mm',1000*clearance);
        statusLabel='No contact';statusColor=[0 .35 .15];
    else
        assert(clearance<0,'The untightened inset must show actual obstacle penetration.');
        distanceLabel=sprintf('Penetration: %.2f mm',-1000*clearance);
        statusLabel='Obstacle penetration';statusColor=[.55 .08 .05];
    end
    % Equal spatial aspect ratio keeps the circular obstacle geometrically true.
    az=axes(fig,'Position',[.775 .345 .205 .43]);hold(az,'on');obstacles(az,P);
    plot(az,xy(:,1),xy(:,2),'Color',cols(ia,:),'LineWidth',1.6);
    plot(az,[pt(1) boundary(1)],[pt(2) boundary(2)],'k--','LineWidth',1.0);
    plot(az,boundary(1),boundary(2),'kx','MarkerSize',5,'LineWidth',1);
    plot(az,pt(1),pt(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cols(ia,:),'MarkerSize',4);
    axis(az,'equal');xlim(az,center(1)+[-half half]);ylim(az,center(2)+[-half half]);styleAxes(az);set(az,'FontSize',8);
    xticks(az,center(1)+[-.8 0 .8]*half);yticks(az,center(2)+[-.8 0 .8]*half);
    xlabel(az,'x (m)');ylabel(az,'y (m)');title(az,'Obstacle boundary','FontWeight','normal','FontSize',9);
    if mode==1,xtickformat(az,'%.2f');ytickformat(az,'%.2f');else,xtickformat(az,'%.3f');ytickformat(az,'%.3f');end
    rectangle(ax,'Position',[center-half,2*half,2*half],'EdgeColor','k','LineStyle',':','LineWidth',.85);
    plot(ax,pt(1),pt(2),'ko','MarkerSize',5,'LineWidth',.8,'HandleVisibility','off');
    annotation(fig,'textbox',[.74 .145 .255 .10],'String',{statusLabel,distanceLabel},'EdgeColor','none','FontName','Times New Roman','FontSize',9,'Color',statusColor,'HorizontalAlignment','center');
    zoomDetails(mode).method=label;zoomDetails(mode).agent=ia;zoomDetails(mode).obstacle=io;
    zoomDetails(mode).physicalTime=D.t(it);zoomDetails(mode).signedClearanceMeters=clearance;
    zoomDetails(mode).trajectoryPoint=pt;zoomDetails(mode).boundaryPoint=boundary;
    zoomDetails(mode).file=[fileName '.fig'];
    saveEditable(fig,outDir,fileName);
end
penetration=-zoomDetails(2).signedClearanceMeters;mn=min(U.obsAll(:));

% Figures 3--5: one inset each, restricted to the untightened violation.
fields={'hInt','hObs','hTree'};
names={'inter_agent_avoidance','static_obstacle_avoidance','connectivity_preservation'};
labels={'$\min_{i<j} h_{ij}^{\mathrm{int}}(y)$','$\min_{i,o} h_{i,o}^{\mathrm{obs}}(y)$','$h^{\mathrm{c}}(y)$'};
mins=zeros(3,2);
for j=1:3
    yt=T.(fields{j});yu=U.(fields{j});[low,ii]=min(yu);tc=U.t(ii);mins(j,:)=[min(yt),low];
    fig=newFigure(names{j},4.65,3.3);ax=axes(fig,'Position',[.14 .16 .82 .75]);hold(ax,'on');
    plot(ax,T.t,yt,'-','Color',[.82 .16 .12],'LineWidth',1.6,'DisplayName','Tightened');
    plot(ax,U.t,yu,'--','Color',[.08 .34 .72],'LineWidth',1.5,'DisplayName','Untightened');
    yline(ax,0,'k-','LineWidth',.85,'HandleVisibility','off');
    xlim(ax,[0 400]);span=max([yt;yu])-min([0;yu]);ylim(ax,[min(0,low)-.07*span,max([yt;yu])+.12*span]);
    styleAxes(ax);xlabel(ax,'Time t (s)');ylabel(ax,labels{j},'Interpreter','latex');
    legend(ax,'Location','northoutside','Orientation','horizontal','Box','off','FontSize',10);
    win=[tc-6 tc+6];small=abs(low);
    insetPositions=[.59 .26 .345 .27;.535 .44 .395 .33;.54 .27 .39 .23];
    az=axes(fig,'Position',insetPositions(j,:));hold(az,'on');
    patch(az,[win(1) win(2) win(2) win(1)],[-1.65*small -1.65*small 0 0],[1 .92 .92],'EdgeColor','none');
    plot(az,U.t,yu,'--','Color',[.08 .34 .72],'LineWidth',1.4);
    yline(az,0,'k-','LineWidth',.9);plot(az,tc,low,'o','Color',[.08 .34 .72],'MarkerFaceColor',[.08 .34 .72],'MarkerSize',3.5);
    xlim(az,win);ylim(az,[-1.65*small,1.0*small]);styleAxes(az);set(az,'FontSize',8.5);
    title(az,'Untightened (\times10^{-3})','FontSize',9,'FontWeight','normal');
    xlabel(az,'t (s)','FontSize',9);az.YAxis.Exponent=0;
    yticklabels(az,compose('%g',az.YTick*1e3));
    saveEditable(fig,outDir,names{j});
end

% Figure 6: position-component convergence, one row and two columns.
% Use physical time from the same tightened Experiment A run.
positionWindow=[0 min(450,T.t(end))];
fig=newFigure('Position-component convergence',7.16,3.35);
positions={[.09 .22 .36 .73],[.60 .22 .36 .73]};
for component=1:2
    ax=subplot(1,2,component,'Parent',fig);set(ax,'Position',positions{component});hold(ax,'on');
    curves=gobjects(5,1);
    for i=1:5
        curves(i)=plot(ax,T.t,T.Y(:,i,component),'Color',cols(i,:),'LineStyle',styles{i},...
            'LineWidth',1.2,'DisplayName',sprintf('Agent %d',i));
    end
    xlim(ax,positionWindow);xticks(ax,0:100:positionWindow(2));
    box(ax,'on');grid(ax,'off');ax.TickDir='in';ax.LineWidth=.7;
    xlabel(ax,'Time t (s)');
    if component==1
        ylabel(ax,'Horizontal Position (m)');ylim(ax,[-3.2 3.2]);yticks(ax,-3:1:3);
    else
        ylabel(ax,'Vertical Position (m)');ylim(ax,[-2.2 2]);yticks(ax,-2:1:2);
    end
    lg=legend(ax,curves,'Location','south','NumColumns',2,'Box','off','FontSize',8.5);
    lg.ItemTokenSize=[12 8];
end
annotation(fig,'textbox',[.2 .012 .6 .07],'String','Position-component convergence',...
    'EdgeColor','none','HorizontalAlignment','center','FontName','Times New Roman','FontSize',10);
saveEditable(fig,outDir,'HB_flow_hebing');

% Figure 7: one error curve only, read from the same full trajectory as above.
E=entryAndError(P,R.C,R.proposed);
tailStart=E.summary.tailStart;tailEnd=E.summary.tailEnd;
select=E.t>=tailStart & E.t<=tailEnd;
assert(E.entryChecks(find(select,1)),'The plotted tail must start inside the certified initialization set.');
fig=newFigure('Local error from the same regional run',4.65,3.3);
ax=axes(fig,'Position',[.16 .17 .80 .77]);
semilogy(ax,E.t(select)-tailStart,E.error(select),'-','Color',[.08 .34 .72],'LineWidth',1.7);
styleAxes(ax);ax.YMinorGrid='off';xlim(ax,[0 tailEnd-tailStart]);ylim(ax,[.65*min(E.error(select)),1.5*max(E.error(select))]);
xlabel(ax,['Elapsed time \tau = t - ' num2str(tailStart) ' (s)']);ylabel(ax,'$\|e(t)\|_2$','Interpreter','latex');
saveEditable(fig,outDir,'exp_2');

report=struct('singleRegionalRun',true,'independentLocalExperiment',false,'controllerSwitch',false,'stateReset',false,...
    'figureCount',7,'insetsPerFigure',[1 1 1 1 1 0 0],'positionComponentTimeWindow',positionWindow,'trajectoryZooms',zoomDetails,'safetyMinimaTightenedUntightened',mins,...
    'obstaclePenetrationMeters',penetration,'obstacleMinimum',mn,'entry',E.summary,...
    'corollary1Claimed',false,'regionalUniformCertificate',false);
fid=fopen(fullfile(outDir,'figure_data_audit.json'),'w');fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));fclose(fid);
writetable(table(E.t(select),E.t(select)-tailStart,E.error(select),E.slow(select),E.fast(select),...
    'VariableNames',{'physical_time','elapsed_time','full_error','slow_error','fast_error'}),fullfile(outDir,'same_run_error_tail.csv'));
writetable(table(T.t,T.hInt,T.hObs,T.hTree,T.lambda2,'VariableNames',{'time','pair_margin','obstacle_margin','connectivity_margin','lambda2'}),fullfile(outDir,'tightened_margins.csv'));
writetable(table(U.t,U.hInt,U.hObs,U.hTree,U.lambda2,'VariableNames',{'time','pair_margin','obstacle_margin','connectivity_margin','lambda2'}),fullfile(outDir,'untightened_margins.csv'));
save(fullfile(outDir,'plot_data.mat'),'T','U','E','-v7');
disp(report.entry);
end

function E=entryAndError(P,C,S)
nt=numel(S.t);slow=zeros(nt,1);fast=slow;energy=slow;dualMax=slow;
Xstar=repmat(reshape(P.ystar.',1,10),5,1);
for j=1:nt
    x=reshape(S.z(j,1:50),10,5).'-Xstar;v=reshape(S.z(j,51:60),2,5).';u=S.z(j,61:70).';aux=S.z(j,71:95).';la=S.z(j,96:140).';
    g=zeros(5,2);y=g;
    for i=1:5,g(i,:)=P.B(i,:)*reshape(x(i,:),2,5).';y(i,:)=x(i,2*i-1:2*i);end
    ef=[u+reshape(g.',10,1);aux];slow(j)=norm([x(:);v(:)]);fast(j)=norm([ef;la]);dualMax(j)=max(abs(la));
    vn=.5*sum(v.^2,'all')+P.k*sum(y.*v,'all')+P.k*P.gamma/2*sum(y.^2,'all')+P.k/2*(sum(x.^2,'all')-sum(y.^2,'all'));
    energy(j)=vn+C.sigma/2*(sum(ef.^2)+C.kappa*sum(la.^2));
end
r1=sqrt(C.rhoComposite/(2*C.mBar));rf=min(C.rFast,sqrt(C.rhoComposite/(C.sigma*max(1,C.kappa))));
checks=slow<r1 & fast<rf & energy<C.rhoComposite & dualMax==0;first=find(checks,1);assert(~isempty(first)&&all(checks(first:end)));
startTime=max(P.tailStart,25*ceil((S.t(first)+10)/25));
start=find(S.t>=startTime,1);finish=find(S.t>=startTime+P.tailDuration,1);
assert(~isempty(finish),'Extend the simulation horizon to include the full local tail.');
error=hypot(slow,fast);tail=start:finish;fit=polyfit(S.t(tail)-S.t(start),log(error(tail)),1);
summary=struct('firstAdmissibleSavedTime',S.t(first),'tailStart',S.t(start),'tailEnd',S.t(finish),...
    'theta1Radius',r1,'fastInitializationRadius',rf,'slowNormAtTailStart',slow(start),'fastNormAtTailStart',fast(start),...
    'dualMaxAtTailStart',dualMax(start),'initialError',error(start),'finalError',error(finish),...
    'relativeFinalError',error(finish)/error(start),'fittedDecayRate',-fit(1),'certifiedNormRate',C.rateBound,...
    'allLaterSavedSamplesAdmissible',true);
E=struct('t',S.t,'slow',slow,'fast',fast,'error',error,'energy',energy,'entryChecks',checks,'summary',summary);
end

function D=geometry(P,S)
t=S.t;nt=numel(t);Y=zeros(nt,5,2);
for i=1:5,j=(i-1)*10+2*i-1;Y(:,i,:)=S.z(:,j:j+1);end
hp=zeros(nt,10);ho=zeros(nt,5,3);hc=zeros(nt,1);la=hc;
for j=1:10,ij=P.edges(j,:);hp(:,j)=sum((Y(:,ij(1),:)-Y(:,ij(2),:)).^2,3)-P.dsafe^2;end
for i=1:5,for j=1:3,ho(:,i,j)=sum((Y(:,i,:)-reshape(P.obs(j,:),1,1,2)).^2,3)-P.rObs(j)^2;end,end
for k=1:nt
    y=squeeze(Y(k,:,:));ds=sum((reshape(y,5,1,2)-reshape(y,1,5,2)).^2,3);W=zeros(5);ok=ds<P.dc^2;
    W(ok)=exp(-ds(ok)./(P.dc^2-ds(ok)));W(1:6:end)=0;L=diag(sum(W,2))-W;ev=sort(eig(L));
    la(k)=ev(2);hc(k)=5*det(L(1:4,1:4))/P.treeScale-P.lambdaLower;
end
D=struct('t',t,'Y',Y,'pairAll',hp,'obsAll',ho,'hInt',min(hp,[],2),'hObs',min(reshape(ho,nt,15),[],2),'hTree',hc,'lambda2',la);
end

function obstacles(ax,P)
a=linspace(0,2*pi,2000);
for j=1:3,fill(ax,P.obs(j,1)+P.rObs(j)*cos(a),P.obs(j,2)+P.rObs(j)*sin(a),[.84 .84 .84],'EdgeColor',[.25 .25 .25],'LineWidth',.8,'HandleVisibility','off');end
end

function fig=newFigure(name,w,h)
fig=figure('Color','w','Units','inches','Position',[1 1 w h],'Name',name,'NumberTitle','off','Renderer','painters');
end

function styleAxes(ax)
grid(ax,'on');box(ax,'on');ax.GridAlpha=.12;ax.LineWidth=.7;ax.TickDir='out';ax.Layer='top';
end

function saveEditable(fig,outDir,name)
set(fig,'Visible','on','WindowStyle','normal','MenuBar','figure','ToolBar','figure','HandleVisibility','on');
set(findall(fig,'Type','axes'),'HitTest','on');drawnow;savefig(fig,fullfile(outDir,[name '.fig']));
set(fig,'Visible','off');exportgraphics(fig,fullfile(outDir,[name '.png']),'Resolution',250);
exportgraphics(fig,fullfile(outDir,[name '.pdf']),'ContentType','vector');
print(fig,fullfile(outDir,[name '.eps']),'-depsc','-painters');close(fig);
end

% ======================================================================
% Embedded source: write_numerical_audit.m
% ======================================================================
function S=write_numerical_audit(R,outDir)
% Readable scalar summary; full time series stay in experiment_data.mat.
if nargin<2,outDir=fullfile(fileparts(mfilename('fullpath')),'results');end
S=struct('parameters',R.P,'certificates',R.C,'sharedController',R.proposed.controller,'corollary1Claimed',false);
if isfield(R,'integrationSettings'),S.integrationSettings=R.integrationSettings;end
drop={'t','ru','distortion','deltaC','Vnom','Vc','error','hInt','hObs','hTree','lambda2'};
if isfield(R,'auditProposed'),S.experimentA=rmfield(R.auditProposed,intersect(fieldnames(R.auditProposed),drop));end
if isfield(R,'auditUntightened'),S.experimentAAblation=rmfield(R.auditUntightened,intersect(fieldnames(R.auditUntightened),drop));end
S.singleRegionalRun=true;S.independentLocalExperiment=false;
S.continuousTimeRegionalCertificate=false;
fid=fopen(fullfile(outDir,'numerical_audit.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(S,'PrettyPrint',true));
end

% ======================================================================
% Embedded source: verify_editable_figures.m
% ======================================================================
function report=verify_editable_figures(folder)
% Open -> change a native line property -> save -> reopen -> verify.
files=dir(fullfile(folder,'*.fig'));report=repmat(struct(),numel(files),1);
checkDir=fullfile(folder,'verification','figure_edit_roundtrip');if ~exist(checkDir,'dir'),mkdir(checkDir);end
for i=1:numel(files)
    file=fullfile(files(i).folder,files(i).name);S=load(file,'-mat');
    root=S.hgS_070000;storedVisible='on';if isfield(root.properties,'Visible'),storedVisible=root.properties.Visible;end
    h=openfig(file,'new','invisible');set(h,'Visible','on','WindowStyle','normal');drawnow;
    assert(strcmp(h.Visible,'on'));assert(strcmp(h.MenuBar,'figure'));assert(strcmp(h.ToolBar,'figure'));
    ln=findall(h,'Type','line');assert(~isempty(ln),'Figure contains no native editable line objects.');
    widths=get(ln,{'LineWidth'});data=get(ln,{'XData','YData'});
    ln(1).Tag='EDITABILITY_ROUNDTRIP';ln(1).LineWidth=widths{1}+.125;plotedit(h,'on');
    target=fullfile(checkDir,files(i).name);savefig(h,target);close(h);
    h2=openfig(target,'new','invisible');edited=findall(h2,'Tag','EDITABILITY_ROUNDTRIP');
    assert(numel(edited)==1 && abs(edited.LineWidth-widths{1}-.125)<1e-12);
    assert(isequaln(edited.XData,data{1,1}) && isequaln(edited.YData,data{1,2}));
    report(i).file=files(i).name;report(i).storedVisible=storedVisible;report(i).nativeLineCount=numel(ln);
    report(i).nativeAxesCount=numel(findall(h2,'Type','axes'));report(i).editSaveReopenPass=true;
    close(h2);assert(strcmp(storedVisible,'on'),'Figure was serialized with Visible=off.');
end
fid=fopen(fullfile(folder,'figure_editability_check.json'),'w');fprintf(fid,'%s',jsonencode(struct('matlabVersion',version,'figures',report,'allPassed',true),'PrettyPrint',true));fclose(fid);
end

