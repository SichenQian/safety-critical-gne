function D=audit_state(P,z,tight)
% Independent physical ECBF assembly; no estimated connectivity rows reused.
X=reshape(z(P.iX),10,5).';V=reshape(z(P.iV),2,5).';U=reshape(z(P.iU),2,5).';
[A,b,nom,aux]=audit_core('rows',P,X,V,tight);Y=aux.Y;
% Recompute physical drift and estimator derivatives directly from edge sums.
W=zeros(5);rho=zeros(5,10);
for i=1:5
 for j=1:5
  s=sum((Y(i,:)-Y(j,:)).^2);
  if i~=j && s<P.dc^2,W(i,j)=exp(-s/(P.dc^2-s));end
  rho(i,:)=rho(i,:)+P.c*W(i,j)*(X(i,:)-X(j,:));
 end
end
L=diag(sum(W,2))-W;drift=-P.gamma*V;Xd=-rho;
for i=1:5,ii=2*i-1:2*i;drift(i,:)=drift(i,:)-rho(i,ii);Xd(i,ii)=V(i,:);end
acc=drift+U;Ac=zeros(P.cConn+5,10);bc=zeros(P.cConn+5,1);q=zeros(P.cConn,1);pairGap=zeros(10,1);weights=pairGap;psiPair=pairGap;
for e=1:10
 i=P.edges(e,1);j=P.edges(e,2);d=Y(i,:)-Y(j,:);v=V(i,:)-V(j,:);h=d*d.'-P.dsafe^2;
 weights(e)=W(i,j);pairGap(e)=norm(d)-P.dsafe;psiPair(e)=2*d*v.'+P.alpha1*h;
 aa=zeros(1,10);aa(2*i-1:2*i)=-2*d;aa(2*j-1:2*j)=2*d;
 bb=-2*d*(drift(i,:)-drift(j,:)).'-2*(v*v.')-2*P.K2*d*v.'-P.K1*h;
 Ac(e,:)=W(i,j)*aa;bc(e)=W(i,j)*(bb+tight*parameter_at(P.DeltaInt,e));q(e)=aa*z(P.iU)+bb;
end
obsGap=zeros(5*P.nObs,1);psiObs=obsGap;
for i=1:5
 for o=1:P.nObs
  r=(i-1)*P.nObs+o;d=Y(i,:)-P.obs(o,:);h=d*d.'-P.rObs(o)^2;
  obsGap(r)=norm(d)-P.rObs(o);psiObs(r)=2*d*V(i,:).'+P.alpha1*h;
  Ac(10+r,2*i-1:2*i)=-2*d;
  bb=-2*d*drift(i,:).'-2*sum(V(i,:).^2)-2*P.K2*d*V(i,:).'-P.K1*h;
  bc(10+r)=bb+tight*parameter_at(P.DeltaObs,r);q(10+r)=Ac(10+r,:)*z(P.iU)+bb;
 end
end
[hc,gc,Hdir]=audit_core('tree',P,Y,V);
hcDirect=(P.N*det(L(1:4,1:4))-512*P.lambdaLower)/512;
Ac(P.cConn,:)=-reshape(gc.',1,10);bb=-sum(gc.*drift,'all')-Hdir-P.K2*sum(gc.*V,'all')-P.K1*hc;
bc(P.cConn)=bb+tight*P.DeltaC;q(P.cConn)=Ac(P.cConn,:)*z(P.iU)+bb;
for i=1:5
 s=sum(V(i,:).^2);tau=max(0,min(1,(4*s-P.nu^2)/(3*P.nu^2)));chi=10*tau^3-15*tau^4+6*tau^5;
 ii=2*i-1:2*i;Ac(P.cConn+i,ii)=chi*V(i,:);bc(P.cConn+i)=-chi*(V(i,:)*nom(ii)+P.theta*s)-(1-chi)*P.deltaP;
end
trueConn=q(P.cConn);hatConn=sum(A(P.rConn,1:10)*z(P.iU)+b(P.rConn))/P.connScale-tight*P.DeltaC;
da=Ac(P.cConn,:)-sum(A(P.rConn,1:10),1)/P.connScale;db=trueConn-hatConn-da*z(P.iU);
eX=X-repmat(reshape(P.ystar.',1,10),5,1);ey=X-repmat(reshape(Y.',1,10),5,1);
ef=[z(P.iU)-nom;z(P.iZ);z(P.iL)];
eigL=sort(eig(L));active=vecnorm(Y(P.edges(:,1),:)-Y(P.edges(:,2),:),2,2)<P.dc;
positive=W(W>0);if isempty(positive),wp=NaN;else,wp=min(positive);end
% Connectivity components use exact nonzero support, independently of eig tol.
reach=logical(eye(5));for ee=find(active).',ii=P.edges(ee,1);jj=P.edges(ee,2);reach(ii,jj)=true;reach(jj,ii)=true;end;for k=1:5,reach=reach|(reach(:,k)&reach(k,:));end
components=0;seen=false(1,5);for i=1:5,if ~seen(i),components=components+1;seen=seen|reach(i,:);end,end
pairRowErr=0;for e=1:10,rr=2*e-1:2*e;pairRowErr=max(pairRowErr,abs(sum(A(rr,:)*z([P.iU P.iZ])+b(rr))-P.barrierScale*(Ac(e,:)*z(P.iU)+bc(e))));end
pairq=q(1:10);pairq(~active)=-Inf;
D=struct('A',A,'b',b,'nom',nom,'Ac',Ac,'bc',bc,'q',q,'W',W,'Y',Y,'V',V,'weights',weights,...
 'pairGap',pairGap,'obsGap',obsGap,'hc',hc,'lambda2',eigL(2),'components',components,'minPositiveWeight',wp,...
 'psiPair',psiPair,'psiObs',psiObs,'psiC',sum(gc.*V,'all')+P.alpha1*hc,...
 'deltaActual',abs(trueConn-hatConn),'deltaBall',abs(db)+P.auditInputRadius*norm(da),...
 'estimateError',norm(ey,'fro'),'slowError',norm([eX(:);V(:)]),'fastLocalError',norm(ef),...
 'trueQmax',[max(pairq),max(q(P.cObs)),q(P.cConn)],'allPairQmax',max(q(1:10)),...
 'flowAgreement',norm([drift(:)-aux.f0(:);Xd(:)-aux.Xdot(:)],inf),...
 'certificateAgreement',abs(hc-hcDirect),'pairSumAgreement',pairRowErr,...
 'dynamicLocalResidual',max(A*z([P.iU P.iZ])+b),'activePairs',nnz(active));
end
