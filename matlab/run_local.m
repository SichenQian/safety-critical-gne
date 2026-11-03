function run_local(name,refine)
if nargin<1,name='local';end
if nargin<2,refine=false;end
root=paper_data_root();folder=fullfile(root,'results',name);if ~exist(folder,'dir'),mkdir(folder);end
F=load(fullfile(root,'results','selected','trajectory.mat'),'P');P=F.P;C=theorem_certificates(P);assert(C.analyticLocalPass);
P.RelTol=2e-10;P.AbsTol=1e-14;P.maxStep=.1;P.T=200;P.outputDt=.1;P.denseDt=.025;
if refine,P.RelTol=2e-11;P.AbsTol=1e-15;P.maxStep=.05;end
P.outputTimes=unique([0;logspace(-9,-2,200).';(.02:.02:1).';(1.1:.1:P.T).']);
% Deterministic non-consensus perturbation of all slow estimates/velocities.
d=sin((1:60)'*sqrt(2));d=d/norm(d)*(.75*C.rTheta1);
zstar=[reshape(repmat(reshape(P.ystar.',1,10),5,1).',50,1);zeros(P.dim-50,1)];
z0=zstar;z0(1:60)=z0(1:60)+d;
X=reshape(z0(P.iX),10,5).';V=reshape(z0(P.iV),2,5).';[A,b,un]=audit_core('rows',P,X,V,true);
z0(P.iU)=un;z0(P.iZ(C.retainedZ(1)))=.5*C.rFast;z0(P.iL)=0;
P.initialState=z0;
eTheta=norm(z0(1:60)-zstar(1:60));ef=norm([z0(P.iU)-un;z0(P.iZ(C.retainedZ));z0(P.iL(C.retainedRows))]);
slow=xslow(P,z0,zstar,C.k);vf=.5*norm([z0(P.iU)-un;z0(P.iZ(C.retainedZ))])^2;
initial=struct('eTheta',eTheta,'rTheta1',C.rTheta1,'eFast',ef,'rFast',C.rFast,'Vc',slow+C.sigma*vf,'rho',C.rho,...
 'maxNominalRow',max(A*[un;zeros(25,1)]+b),'maxActualRow',max(A*z0([P.iU P.iZ])+b),...
 'allDualZero',all(z0(P.iL)==0),'omittedPairs',C.omittedPairs,'certified',eTheta<C.rTheta1 && ef<C.rFast && slow+C.sigma*vf<C.rho);
gn=A*[un;zeros(25,1)]+b;ga=A*z0([P.iU P.iZ])+b;
initial.maxRetainedNominalRow=max(gn(C.retainedRows));initial.maxRetainedActualRow=max(ga(C.retainedRows));
assert(initial.certified && initial.maxRetainedNominalRow<0 && initial.maxRetainedActualRow<0);write_json(fullfile(folder,'certificate.json'),C);write_json(fullfile(folder,'initial.json'),initial);
S=simulate_local_error(P,zstar);save(fullfile(folder,'trajectory.mat'),'P','S','C','initial','zstar','-v7');
audit_final(name,'full');validate_local_certificate(name);
end
function v=xslow(P,z,zs,k)
X=reshape(z(P.iX)-zs(P.iX),10,5).';vv=z(P.iV);yy=zeros(10,1);
for i=1:5,yy(2*i-1:2*i)=X(i,2*i-1:2*i).';X(i,2*i-1:2*i)=0;end
v=.5*(vv.'*vv)+k*(yy.'*vv)+k*P.gamma/2*(yy.'*yy)+k/2*sum(X.^2,'all');
end
