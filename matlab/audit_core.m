function varargout = audit_core(action,varargin)
%AUDIT_CORE Snapshot-derived equations. No model/parameter changes by default.
% See CHANGELOG.md for instrumentation and the explicit stress-only overrides.
switch action
case 'parameters'
 P=parameters();P.k=sqrt(P.c*(P.N-1)+1.15^2+sum(P.B(1,:).^2)/(2*(P.c*P.lambdaProof/2-norm(P.B(1,:))-(1.15+norm(P.B(1,:)))^2/2)));
 P.solver='ode15s';P.captureMesh=true;P.denseDt=.025;P.storeSolutions=true;P.maxWallSeconds=900;
 P.estimateOffset=zeros(5,10);varargout={P};
case 'configure',varargout={configure(varargin{:})};
case 'certificates',varargout={theorem_certificates(varargin{:})};
case 'simulate',varargout={simulate(varargin{:})};
case 'rows',[varargout{1:nargout}]=rows(varargin{:});
case 'graph',[varargout{1:nargout}]=graph(varargin{:});
case 'tree',[varargout{1:nargout}]=treeBarrier(varargin{:});
case 'det',[varargout{1:nargout}]=detDerivatives(varargin{:});
case 'qp',[varargout{1:nargout}]=solveScaledQP(varargin{:});
case 'verify',varargout={verifyImplementation(varargin{:})};
otherwise,error('Unknown action');
end
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
P.treeScale=512;P=configure(P);
end

function P=configure(P)
P.nObs=size(P.obs,1);assert(numel(P.rObs)==P.nObs);
P.nLambda=30+5*P.nObs;P.dim=95+P.nLambda;P.iL=96:P.dim;
P.rPair=1:20;P.rObsRows=21:20+5*P.nObs;
P.rConn=21+5*P.nObs:25+5*P.nObs;P.rDiss=26+5*P.nObs:P.nLambda;
P.cObs=11:10+5*P.nObs;P.cConn=11+5*P.nObs;P.cDiss=P.cConn+(1:5);
P.alpha1=(P.K2-sqrt(P.K2^2-4*P.K1))/2;P.alpha2=P.K2-P.alpha1;
P.lambdaProof=P.lambdaLower;
lt=norm(P.B(1,:));mu=min(eig(P.B));ll=norm(P.B,2);ky=P.c*P.lambdaLower/2-lt-(ll+lt)^2/(2*mu);
if ky>0,P.k=sqrt(P.c*(P.N-1)+ll^2/mu+lt^2/(2*ky));else,P.k=NaN;end
end

function res=simulate(P,tight)
% Integrate each smooth face and locate changes of the projected dual flow
% with events. This avoids finite-difference chattering at lambda=0.
X0=repmat(reshape(P.y0.',1,10),5,1);V0=P.v0;X0=X0+P.estimateOffset;
if isfield(P,'initialState'),X0=reshape(P.initialState(P.iX),10,5).';V0=reshape(P.initialState(P.iV),2,5).';end
[A,b,nom]=rows(P,X0,V0,tight);H=diag([ones(10,1);P.xi*ones(25,1)]);
opt=optimoptions('quadprog','Display','off','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
if ~isfield(P,'initialState')
[w,exitflag,lm]=solveScaledQP(P,A,b,nom);
if isfield(P,'exactInitialization') && P.exactInitialization && max(A*[nom;zeros(25,1)]+b)<0
 w=[nom;zeros(25,1)];exitflag=1;lm.ineqlin=zeros(P.nLambda,1);
end
assert(exitflag>0,'Initial QP solver did not confirm a solution; independent feasibility check required.');la=max(0,lm.ineqlin);la(A*w+b<-1e-7)=0;
else,w=P.initialState([P.iU P.iZ]);la=P.initialState(P.iL);exitflag=1;end
zz0=[reshape(X0.',50,1);reshape(V0.',10,1);w;la];active=la>1e-9;zz0(P.iL(~active))=0;
if isfield(P,'initialState')
 zz0=P.initialState;la=zz0(P.iL);assert(all(la>=0));active=la>0;
end
gridT=(0:P.outputDt:P.T).';if isfield(P,'outputTimes'),gridT=P.outputTimes(:);end
allZ=nan(numel(gridT),P.dim);allZ(1,:)=zz0.';solutions={};
t0=0;start=tic;lastPrint=-100;lastRhs=0;nEvents=0;eventLog=zeros(0,2);meshT=[];meshZ=[];denseT=[];denseZ=[];crossings=zeros(0,7);rhsCalls=0;jacCalls=0;stepCount=0;
fprintf('Integrating projected flow by smooth faces: tight=%d\n',tight);
fprintf('Initial filter distortion %.9g, active multipliers %d\n',norm(w(1:10)-nom),nnz(active));
while t0<P.T
    opts=odeset('RelTol',P.RelTol,'AbsTol',P.AbsTol,'MaxStep',P.maxStep,...
        'Jacobian',@jacobian,'Events',@events,'OutputFcn',@progress);
    solver=str2func(P.solver);sol=solver(@rhs,[t0 P.T],zz0,opts);tend=sol.x(end);
    stepCount=stepCount+numel(sol.x)-1;
    % deval uses x/y/idata; retain extdata as an empty marker for modern
    % interpolation format. Strip nested RHS handles before serialization
    % to avoid recursively capturing every earlier segment/workspace.
    if P.storeSolutions,cleanSol=sol;cleanSol.extdata=struct();solutions{end+1}=cleanSol;end %#ok<AGROW>
    if P.captureMesh
      meshT=[meshT;sol.x(:)];meshZ=[meshZ;sol.y.']; %#ok<AGROW>
      % Dense physical grid plus logarithmically close face-event samples.
      % All internal accepted states remain available separately in meshZ.
      tt=(ceil(t0/P.denseDt)*P.denseDt:P.denseDt:tend).';off=[0,1e-8,1e-6,1e-4,.001,.01];near=[t0+off;tend-off];
      tt=unique([tt;near(:);t0;tend]);tt=tt(tt>=t0 & tt<=tend);
      denseT=[denseT;tt];denseZ=[denseZ;deval(sol,tt).']; %#ok<AGROW>
      for ee=1:P.m
        ij=P.edges(ee,:);ci=(ij(1)-1)*10+2*ij(1)-1; cj=(ij(2)-1)*10+2*ij(2)-1;
        gap=vecnorm(sol.y(ci:ci+1,:)-sol.y(cj:cj+1,:),2,1)-P.dc;
        kk=find(gap(1:end-1).*gap(2:end)<0);
        for crossIndex=kk
          tr=fzero(@(t)crossgap(t,sol,ci,cj),[sol.x(crossIndex),sol.x(crossIndex+1)]);zr=deval(sol,tr);
          dv=reshape(zr(P.iV),2,5);dr=zr(ci:ci+1)-zr(cj:cj+1);vr=dv(:,ij(1))-dv(:,ij(2));
          radial=dr.'*vr/P.dc;psi=2*dr.'*vr+P.alpha1*(P.dc^2-P.dsafe^2);
          crossings=[crossings;tr,ij,radial,norm(vr),psi,double(radial<0)]; %#ok<AGROW>
        end
      end
    end
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
res.internalT=meshT;res.internalZ=meshZ;res.denseT=denseT;res.denseZ=denseZ;res.sensingCrossings=crossings;res.rhsCalls=rhsCalls;res.jacCalls=jacCalls;res.internalSteps=stepCount;res.solver=P.solver;res.controller=controller(P);if ~tight,res.controller.DeltaInt=0;res.controller.DeltaObs=0;res.controller.DeltaC=0;end
res.solutions=solutions;
fprintf('Finished tight=%d, t=%.0f, events=%d, wall=%.1f s\n',tight,P.T,nEvents,toc(start));
    function dz=rhs(tt,zz)
        rhsCalls=rhsCalls+1;if toc(start)>P.maxWallSeconds,error('audit:Timeout','Integration exceeded declared wall time cap.');end
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
        jacCalls=jacCalls+1;
        J=zeros(P.dim);X=reshape(zz(P.iX),10,5).';V=reshape(zz(P.iV),2,5).';[AA,~]=rows(P,X,V,tight);
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
        direction=ones(P.nLambda,1);direction(active)=-1;isterminal=ones(P.nLambda,1);
        if tt<=t0+4*eps(max(1,t0))
            value(active)=max(value(active),1e-12);value(~active)=min(value(~active),-1e-12);
        end
        % An identically zero inactive collision row has projected derivative
        % zero, not an activation event. Prevent spurious zero-root toggles.
        zeroInactive=false(P.nLambda,1);
        zeroInactive(1:20)=~active(1:20) & all(AA(1:20,:)==0,2) & bb(1:20)==0;
        value(zeroInactive)=-1;
    end
    function g=crossgap(t,ss,ci,cj)
        q=deval(ss,t);g=norm(q(ci:ci+1)-q(cj:cj+1))-P.dc;
    end
    function stop=progress(tt,~,flag)
        stop=false;if isempty(flag)&&~isempty(tt)&&tt(end)-lastPrint>=100
            fprintf('  tight=%d t=%.3f / %.0f, wall=%.1f s\n',tight,tt(end),P.T,toc(start));lastPrint=tt(end);
        end
    end
end

function [A,b,nom,aux]=rows(P,X,V,tight)
Y=zeros(5,2);for i=1:5,Y(i,:)=X(i,2*i-1:2*i);end
centered=isfield(P,'errorCoordinates') && P.errorCoordinates;
if centered,Y=Y+P.ystar;end
[W,L]=graph(P,Y);rho=P.c*L*X;Xdot=-rho;
f0=-P.gamma*V;grad=zeros(5,2);
for i=1:5
    ii=2*i-1:2*i;Xdot(i,ii)=V(i,:);f0(i,:)=f0(i,:)-rho(i,ii);
    localY=reshape(X(i,:),2,5).';grad(i,:)=P.B(i,:)*localY;
    if ~centered,grad(i,:)=grad(i,:)-P.p(i,:);end
end
nom=-reshape(grad.',10,1);A=zeros(P.nLambda,35);b=zeros(P.nLambda,1);hInt=zeros(10,1);hObs=zeros(5*P.nObs,1);rr=0;
for e=1:10
    i=P.edges(e,1);j=P.edges(e,2);d=Y(i,:)-Y(j,:);v=V(i,:)-V(j,:);h=sum(d.^2)-P.dsafe^2;hInt(e)=h;
    for side=1:2
        if side==1,a=i;sgn=1;other=2*e;own=2*e-1;else,a=j;sgn=-1;other=2*e-1;own=2*e;end
        g=sgn*2*d;rr=rr+1;A(rr,2*a-1:2*a)=-W(i,j)*g;
        A(rr,10+own)=W(i,j);A(rr,10+other)=-W(i,j);
        b(rr)=W(i,j)*(-g*f0(a,:).'-P.K2*g*V(a,:).'-sum(v.^2)-.5*P.K1*h+.5*tight*parameter_at(P.DeltaInt,e));
    end
end
for i=1:5
    for o=1:P.nObs
        d=Y(i,:)-P.obs(o,:);h=sum(d.^2)-P.rObs(o)^2;hObs((i-1)*P.nObs+o)=h;g=2*d;
        rr=rr+1;A(rr,2*i-1:2*i)=-g;
        b(rr)=-g*f0(i,:).'-2*sum(V(i,:).^2)-P.K2*g*V(i,:).'-P.K1*h+tight*parameter_at(P.DeltaObs,(i-1)*P.nObs+o);
    end
end
for i=1:5
    Yi=reshape(X(i,:),2,5).';if centered,Yi=Yi+P.ystar;end
    Vi=reshape(Xdot(i,:),2,5).';
    [hc,gc,Hdir]=treeBarrier(P,Yi,Vi);
    rr=rr+1;A(rr,2*i-1:2*i)=-gc(i,:);A(rr,31:35)=L(i,:);
    b(rr)=-gc(i,:)*f0(i,:).'-P.K2*gc(i,:)*V(i,:).'-(Hdir+P.K1*hc)/5+tight*P.DeltaC/5;
end
% Positive barrier normalization: pair/obstacle use barrierScale, connectivity
% uses connScale. Each Delta has the same scaling as its barrier.
% Mismatch matrices retain the manuscript form in both experiments/ablations.
A(1:20+5*P.nObs,1:10)=P.barrierScale*A(1:20+5*P.nObs,1:10);b(1:20+5*P.nObs)=P.barrierScale*b(1:20+5*P.nObs);
A(P.rConn,1:10)=P.connScale*A(P.rConn,1:10);b(P.rConn)=P.connScale*b(P.rConn);
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

function [w,flag,lm]=solveScaledQP(P,A,b,nom)
% Exact variable/constraint change for the OFFLINE QP only.
% It leaves the original objective, feasible set, solution and KKT multiplier
% mapping unchanged; the online primal-dual differential equation is untouched.
beta=sqrt(P.barrierScale);bc=sqrt(P.connScale);
T=[ones(10,1);beta*ones(20,1);bc*ones(5,1)];row=[ones(20+5*P.nObs,1)/beta;ones(5,1)/bc;ones(5,1)];
At=(A.*T.').*row;bt=b.*row;Ht=diag([ones(10,1);P.xi*beta^2*ones(20,1);P.xi*bc^2*ones(5,1)]);
opt=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex',...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000);
[x,~,flag,~,ll]=quadprog(Ht,[-nom;zeros(25,1)],At,-bt,[],[],[],[],[nom;zeros(25,1)],opt);
w=T.*x;lm=struct('ineqlin',row.*ll.ineqlin);
end

function S=controller(P)
fields={'c','gamma','epsilon','xi','K1','K2','DeltaInt','DeltaObs','DeltaC','nu','theta','deltaP','dc','dsafe','lambdaLower','barrierScale','connScale'};
S=struct();for j=1:numel(fields),S.(fields{j})=P.(fields{j});end
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
