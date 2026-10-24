function D = diagnose_tightening_margins(source, outDir)
%DIAGNOSE_TIGHTENING_MARGINS Theorem 2 trajectory-based numerical diagnostics.
%   D = DIAGNOSE_TIGHTENING_MARGINS reads results/experiment_data.mat.
%   D = DIAGNOSE_TIGHTENING_MARGINS(MATFILE,OUTDIR) selects saved data/output.
%   D = DIAGNOSE_TIGHTENING_MARGINS(RESULT,OUTDIR) reuses a fresh result and
%   its already solved regularized QPs. Older MAT files lacking QP optimizers
%   are handled with the SAME rows and solver as reproduce_paper, never with
%   the nominal input. No integration, plotting, or source-data writes occur.
%   Default output: reproduced_results/tightening_diagnostics/.
%
%   All reported coefficients and margins use physical (unscaled) row units.
%   The connectivity barrier already includes the original division by 512.
%   delta_c is maximized over ||u|| <= P.auditInputRadius at each saved state;
%   maximizing these values over samples is NOT a maximum over Theta_R.

repo=fileparts(mfilename('fullpath'));
if nargin<1||isempty(source),source=fullfile(repo,'results','experiment_data.mat');end
if nargin<2||isempty(outDir),outDir=fullfile(repo,'reproduced_results','tightening_diagnostics');end
if isstruct(source)
    R=source;sourceLabel='in-memory reproduction result';
else
    loaded=load(source,'result');R=loaded.result;sourceLabel=char(source);
end
P=R.P;C=R.C;S=R.proposed;Aold=R.auditProposed;
assert(S.tight && P.N==5 && P.n==2 && P.nU==10 && P.nZ==25, ...
    'This diagnostic requires the existing five-agent tightened example.');
assert(P.barrierScale>0 && P.connScale>0 && P.auditInputRadius>0, ...
    'Positive row scales and the existing input-ball radius are required.');
assert(all(isfinite(S.z),'all'),'The saved trajectory contains nonfinite states.');
F=reproduce_paper([],'diagnostic_helpers');
t=S.t(:);ns=numel(t);H=diag([ones(P.nU,1);P.xi*ones(P.nZ,1)]);
reuse=isfield(Aold,'qpReference') && isequal(Aold.t(:),t);
if reuse
    Q=Aold.qpReference;
    assert(isequal(size(Q.wStar),[ns,P.nU+P.nZ]) && ...
        isequal(size(Q.lambdaStar),[ns,P.nLambda]) && numel(Q.exitflag)==ns, ...
        'Saved QP reference dimensions do not match the trajectory.');
else
    Q=struct('wStar',nan(ns,P.nU+P.nZ),'lambdaStar',nan(ns,P.nLambda), ...
        'exitflag',nan(ns,1));
end
fprintf('Theorem 2 diagnostic: %d saved states; reuse QP optimizers = %d\n',ns,reuse);
tracking=nan(ns,1);deltaActual=tracking;deltaBall=tracking;qpResidual=tracking;
stationarity=tracking;complementarity=tracking;dualViolation=tracking;
qMax=nan(ns,3);coefficientMax=nan(ns,3);inputNorm=tracking;optimizerNorm=tracking;
localResidual=tracking;distortion=tracking;rowAgreement=tracking;
physicalInside=false(ns,1);ballEndpointError=tracking;activePairs=zeros(ns,1);
scale=[P.barrierScale,P.barrierScale,P.connScale];
Delta=[P.DeltaInt,P.DeltaObs,P.DeltaC];
% Reuse the declared geometry and the existing GLOBAL matrix-tree derivative
% bound. No empirical maximum is substituted for these coefficient bounds.
cornersDistance=zeros(size(P.obs,1),1);
for j=1:size(P.obs,1)
    cornersDistance(j)=hypot(max(abs(P.regionBox(1,:)-P.obs(j,1))), ...
        max(abs(P.regionBox(2,:)-P.obs(j,2))));
end
barA=[2*P.regionDiameter,2*max(cornersDistance),C.Gh];
assert(all(isfinite(barA)) && all(barA>0),'Coefficient bounds are unavailable.');
for k=1:ns
    if mod(k-1,2000)==0,fprintf('  tightening diagnostic %d / %d\n',k,ns);end
    z=S.z(k,:).';X=reshape(z(P.iX),10,5).';V=reshape(z(P.iV),2,5).';
    [A,b,nom,aux]=F.rows(P,X,V,true);u=z(P.iU);
    if ~reuse
        [ws,flag,lm]=F.solveQP(P,A,b,nom);
        Q.exitflag(k)=flag;
        if numel(ws)==P.nU+P.nZ && numel(lm.ineqlin)==P.nLambda
            Q.wStar(k,:)=ws.';Q.lambdaStar(k,:)=lm.ineqlin.';
        end
    end
    ws=Q.wStar(k,:).';lambda=Q.lambdaStar(k,:).';
    if Q.exitflag(k)>0 && all(isfinite([ws;lambda]))
        tracking(k)=norm(u-ws(1:P.nU));
        qpResidual(k)=max(A*ws+b);
        stationarity(k)=norm(H*ws-[nom;zeros(P.nZ,1)]+A.'*lambda,inf);
        complementarity(k)=norm(lambda.*(A*ws+b),inf);
        dualViolation(k)=max(0,-min(lambda));
        optimizerNorm(k)=norm(ws(1:P.nU));distortion(k)=norm(ws(1:P.nU)-nom);
    end
    inputNorm(k)=norm(u);localResidual(k)=max(A*z([P.iU,P.iZ])+b);
    acc=aux.f0+reshape(u,2,P.N).';
    qInt=zeros(P.m,1);pairNorm=zeros(P.m,1);pairWeight=zeros(P.m,1);
    pairRowError=0;
    for e=1:P.m
        i=P.edges(e,1);j=P.edges(e,2);d=aux.Y(i,:)-aux.Y(j,:);v=V(i,:)-V(j,:);
        qInt(e)=-2*d*(acc(i,:)-acc(j,:)).'-2*sum(v.^2) ...
            -2*P.K2*d*v.'-P.K1*aux.hInt(e);
        pairNorm(e)=2*norm(d);pairWeight(e)=aux.W(i,j);
        rr=2*e-1:2*e;
        % Compare weighted rows, without division by small/zero graph weights.
        pairRowError=max(pairRowError,abs(sum(A(rr,1:P.nU)*u+b(rr))- ...
            P.barrierScale*pairWeight(e)*(qInt(e)+P.DeltaInt)));
    end
    qObs=zeros(numel(aux.hObs),1);obsNorm=qObs;
    for i=1:P.N
        for j=1:size(P.obs,1)
            r=(i-1)*size(P.obs,1)+j;d=aux.Y(i,:)-P.obs(j,:);
            qObs(r)=-2*d*acc(i,:).'-2*sum(V(i,:).^2) ...
                -2*P.K2*d*V(i,:).'-P.K1*aux.hObs(r);
            obsNorm(r)=2*norm(d);
        end
    end
    qConn=-sum(aux.gc.*acc,'all')-aux.Hdir-P.K2*sum(aux.gc.*V,'all')-P.K1*aux.hC;
    qHat=sum(A(36:40,1:P.nU)*u+b(36:40))/P.connScale-P.DeltaC;
    da=-reshape(aux.gc.',1,P.nU)-sum(A(36:40,1:P.nU),1)/P.connScale;
    db=qConn-qHat-da*u;
    deltaActual(k)=abs(qConn-qHat);
    deltaBall(k)=abs(db)+P.auditInputRadius*norm(da);
    % The affine discrepancy has its exact input-ball supremum at an endpoint.
    if norm(da)>0
        endpoint=P.auditInputRadius*da.'/norm(da);
        endpointMax=max(abs([da*endpoint+db,-da*endpoint+db]));
    else
        endpointMax=abs(db);
    end
    ballEndpointError(k)=abs(endpointMax-deltaBall(k));
    active=pairWeight>0;activePairs(k)=nnz(active);
    if any(active),qMax(k,1)=max(qInt(active));else,qMax(k,1)=-Inf;end
    qMax(k,2:3)=[max(qObs),qConn];
    coefficientMax(k,:)=[max(pairNorm),max(obsNorm), ...
        max(vecnorm(A(36:40,1:P.nU),2,2))/P.connScale];
    rowAgreement(k)=max([pairRowError; ...
        abs(A(21:35,1:P.nU)*u+b(21:35)-P.barrierScale*(qObs+P.DeltaObs)); ...
        abs(sum(A(36:40,P.nU+1:end),1)).']);
    physicalInside(k)=max(pairNorm)/2<=P.regionDiameter && ...
        max(vecnorm(V,2,2))<=P.regionSpeedMax && ...
        all(aux.Y>=P.regionBox(:,1).' & aux.Y<=P.regionBox(:,2).','all');
end
% Verify the half-pair and 1/N connectivity tightening factors in the actual
% row builder at the initial, maximum-tracking, and final saved states.
[~,peak]=max(tracking);tighteningRowError=0;
for k=unique([1,peak,ns])
    z=S.z(k,:).';X=reshape(z(P.iX),10,5).';V=reshape(z(P.iV),2,5).';
    [A,b,~,aux]=F.rows(P,X,V,true);[A0,b0]=F.rows(P,X,V,false);
    weights=arrayfun(@(e)aux.W(P.edges(e,1),P.edges(e,2)),(1:P.m).');
    expected=[P.barrierScale*repelem(weights,2)*P.DeltaInt/2; ...
        P.barrierScale*P.DeltaObs*ones(15,1); ...
        P.connScale*P.DeltaC/P.N*ones(P.N,1);zeros(P.N,1)];
    tighteningRowError=max([tighteningRowError;abs(A(:)-A0(:));abs(b-b0-expected)]);
end
% Numerical tolerances assess solver/algebra residuals, not safety budgets.
tol=1e-8;validQP=all(Q.exitflag>0) && all(isfinite(tracking));
if validQP,ru=max(tracking);else,ru=NaN;end
dc=max(deltaBall);lhs=[sqrt(2)*barA(1)*ru,barA(2)*ru,sqrt(P.N)*barA(3)*ru+dc];
margin=Delta-lhs;
instantMargins=Delta-tracking*[sqrt(2)*barA(1),barA(2),sqrt(P.N)*barA(3)];
instantMargins(:,3)=instantMargins(:,3)-deltaBall;
checks=validQP && all(qpResidual<=tol) && all(stationarity<=tol) && ...
    all(complementarity<=tol) && all(dualViolation<=tol);
referenceTrackingError=NaN;
if isequal(Aold.t(:),t),referenceTrackingError=max(abs(tracking-Aold.ru));end
summary=struct('scope','trajectory-based numerical diagnostics', ...
    'source',sourceLabel,'samples',ns,'timeInterval',[t(1),t(end)], ...
    'sampleSpacing',P.outputDt,'qpOptimizersReused',reuse,'qpSolveCount',ns*double(~reuse), ...
    'regularizationXi',P.xi,'inputBallRadius',P.auditInputRadius, ...
    'coefficientBoundsPhysical',barA,'rowScales',scale, ...
    'coefficientMaxSampled',max(coefficientMax,[],1),'r_u_sampled',ru, ...
    'delta_c_actual_input_sampled',max(deltaActual), ...
    'delta_c_input_ball_sampled_states',dc,'tighteningPhysical',Delta, ...
    'lhsPhysical',lhs,'marginsPhysical',margin,'marginsScaled',scale.*margin, ...
    'nonnegativeMargins',isfinite(margin)&margin>=0, ...
    'allQPExitFlagsPositive',all(Q.exitflag>0),'allQPNumericalChecksPass',checks, ...
    'numericalResidualTolerance',tol,'maxQPResidual',max(qpResidual), ...
    'maxQPStationarityResidual',max(stationarity), ...
    'maxQPComplementarityResidual',max(complementarity), ...
    'maxQPDualViolation',max(dualViolation), ...
    'maxActualECBFResidualPhysical',max(qMax,[],1), ...
    'maxActualECBFResidualScaled',scale.*max(qMax,[],1), ...
    'allActualECBFResidualsNonpositive',all(qMax<=0,'all'), ...
    'maxDynamicLocalRowResidual',max(localResidual), ...
    'maxInputNorm',max(inputNorm),'maxOptimizerInputNorm',max(optimizerNorm), ...
    'maxOptimizerNominalDifference',max(distortion), ...
    'physicalRegionContainmentSampled',all(physicalInside), ...
    'inputBallContainmentSampled',all(inputNorm<=P.auditInputRadius), ...
    'coefficientBoundsContainSamples',all(coefficientMax<=barA+tol,'all'), ...
    'maxRowAssemblyError',max(rowAgreement),'maxTighteningScalingError',tighteningRowError, ...
    'maxInputBallEndpointError',max(ballEndpointError), ...
    'trackingDifferenceFromOriginalAudit',referenceTrackingError, ...
    'regionalUniformCertificate',false,'continuousTimeRegionalCertificate',false);
summary.sampleSpacingRange=[min(diff(t)),max(diff(t))];
if diff(summary.sampleSpacingRange)>1e-9,summary.sampleSpacing=NaN;end
summary.sampledChecksPass=checks && all(summary.nonnegativeMargins) && ...
    summary.allActualECBFResidualsNonpositive && summary.physicalRegionContainmentSampled && ...
    summary.inputBallContainmentSampled && summary.coefficientBoundsContainSamples && ...
    summary.maxRowAssemblyError<=tol && tighteningRowError<=tol && ...
    summary.maxInputBallEndpointError<=tol;
summary.unverifiedConditions={ ...
    'QP feasibility for every theta in the compact regional set Theta_R', ...
    'A continuous-time tracking bound and containment in Theta_R and U_R^u', ...
    'Connectivity discrepancy over all Theta_R times U_R^u', ...
    'The a priori sensing-boundary relative-speed bound in Assumption 2'};
D=struct('summary',summary,'t',t,'trackingErrorNorm',tracking, ...
    'uStar',Q.wStar(:,1:P.nU),'qpReference',Q,'deltaActual',deltaActual, ...
    'deltaInputBall',deltaBall,'instantMarginsPhysical',instantMargins, ...
    'actualECBFResidualMaxPhysical',qMax,'qpResidual',qpResidual, ...
    'stationarityResidual',stationarity,'complementarityResidual',complementarity, ...
    'activeCollisionPairCount',activePairs);
report=evalc('print_summary(summary)');fprintf('%s',report);
if ~isfolder(outDir),mkdir(outDir);end
save(fullfile(outDir,'tightening_diagnostics.mat'),'D','-v7');
write_text(fullfile(outDir,'tightening_summary.txt'),report);
write_text(fullfile(outDir,'tightening_summary.json'),jsonencode(summary,'PrettyPrint',true));
end

function print_summary(S)
fprintf('\nTheorem 2: trajectory-based numerical diagnostics\n');
fprintf('Samples: %d; t = [%g, %g] s; physical/unscaled row units\n',S.samples,S.timeInterval);
fprintf('r_u (sampled)                 = %.17g\n',S.r_u_sampled);
fprintf('delta_c (sampled input ball)  = %.17g   [||u|| <= %g]\n', ...
    S.delta_c_input_ball_sampled_states,S.inputBallRadius);
fprintf('delta_c (actual input only)   = %.17g\n',S.delta_c_actual_input_sampled);
fprintf('%-14s %12s %12s %12s %12s %8s\n','Condition','bar_a','Delta','LHS','margin','>= 0');
labels={'inter-agent','obstacle','connectivity'};
for i=1:3
    fprintf('%-14s %12.8g %12.8g %12.8g %12.8g %8d\n',labels{i}, ...
        S.coefficientBoundsPhysical(i),S.tighteningPhysical(i),S.lhsPhysical(i), ...
        S.marginsPhysical(i),S.nonnegativeMargins(i));
end
fprintf('QP feasibility/KKT numerical checks: %d (tolerance %.1g)\n', ...
    S.allQPNumericalChecksPass,S.numericalResidualTolerance);
fprintf('Max original ECBF residuals (q <= 0): %.12g  %.12g  %.12g\n', ...
    S.maxActualECBFResidualPhysical);
fprintf('All sampled checks pass: %d; uniform regional certificate: 0\n',S.sampledChecksPass);
fprintf('Finite saved samples do not certify Theorem 2 on the full compact region or between samples.\n');
end

function write_text(path,content)
fid=fopen(path,'w');assert(fid>=0,'Cannot write diagnostic report.');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',content);
end
