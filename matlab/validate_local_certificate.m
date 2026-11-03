function validate_local_certificate(name)
if nargin<1,name='local';end
root=paper_data_root();folder=fullfile(root,'results',name);
D=load(fullfile(folder,'trajectory.mat'));P=D.P;S=D.S;C=D.C;zs=D.zstar;
data=zeros(numel(S.t),5);
for k=1:numel(S.t)
 z=S.z(k,:).';X=reshape(z(P.iX),10,5).';V=reshape(z(P.iV),2,5).';
 [A,b,un]=audit_core('rows',P,X,V,true);eX=reshape(z(P.iX)-zs(P.iX),10,5).';y=zeros(10,1);
 for i=1:5,y(2*i-1:2*i)=eX(i,2*i-1:2*i).';eX(i,2*i-1:2*i)=0;end
 v=z(P.iV);ef=[z(P.iU)-un;z(P.iZ(C.retainedZ));z(P.iL(C.retainedRows))];
 vs=.5*(v.'*v)+C.k*(y.'*v)+C.k*P.gamma/2*(y.'*y)+C.k/2*sum(eX.^2,'all');
 vf=.5*norm([z(P.iU)-un;z(P.iZ(C.retainedZ))])^2+C.kappa/2*norm(z(P.iL(C.retainedRows)))^2;
 residual=A*z([P.iU P.iZ])+b;
 data(k,:)=[norm(z(1:60)-zs(1:60)),norm(ef),vs+C.sigma*vf,max(residual(C.retainedRows)),max(abs(z(P.iL)))];
end
report=struct('slowMaximum',max(data(:,1)),'slowDomainRadius',C.rTheta0,'fastMaximum',max(data(:,2)),'fastTubeRadius',C.rLoc,...
 'VcMaximum',max(data(:,3)),'rho',C.rho,'maxRetainedRow',max(data(:,4)),'maxAbsDual',max(data(:,5)),...
 'initialCertified',D.initial.certified,'analyticLocalPass',C.analyticLocalPass,...
 'allSavedStatesInTube',all(data(:,1)<C.rTheta0 & data(:,2)<C.rLoc & data(:,3)<C.rho));
report.outputGridSamples=numel(S.t);report.scope='Regular saved output grid; omitted identically-zero rows are excluded from retained-row slack.';
write_json(fullfile(folder,'tube_validation.json'),report);assert(report.allSavedStatesInTube && report.analyticLocalPass && report.maxRetainedRow<0);
end
