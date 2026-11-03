function C=theorem_certificates(P)
% Explicit conservative instantiation of the CURRENT v2 Appendix.
% Analytic inequalities evaluated in double precision, not interval arithmetic.
% This is a LOCAL certificate, never a uniform regional safety certificate.
N=P.N; C.mu=min(eig(P.B));C.l=norm(P.B,2);
C.tildeL=max(vecnorm(P.B,2,2));
C.cLower=(2*C.tildeL+(C.l+C.tildeL)^2/C.mu)/P.lambdaLower;
C.kappaY=P.c*P.lambdaLower/2-C.tildeL-(C.l+C.tildeL)^2/(2*C.mu);
C.k=sqrt(P.c*(N-1)+C.l^2/C.mu+C.tildeL^2/(2*C.kappaY));
C.gammaLower=2*C.k;
C.nominalCoefficients=[C.k*C.mu/(4*N),C.k*C.kappaY/2,P.gamma-2*C.k];
C.alphaNom=min(C.nominalCoefficients);
Pk=[C.k*P.gamma C.k;C.k 1];
C.mUnder=.5*min([eig(Pk);C.k]);C.mBar=.5*max([eig(Pk);C.k]);
C.theorem1=P.c>C.cLower && P.gamma>C.gammaLower && C.alphaNom>0;
X=repmat(reshape(P.ystar.',1,10),5,1);V=zeros(5,2);
[A,b,un,aux]=audit_core('rows',P,X,V,true);
dist=vecnorm(P.ystar(P.edges(:,1),:)-P.ystar(P.edges(:,2),:),2,2);
keep=dist<P.dc;C.omittedPairs=find(~keep);C.retainedRows=[reshape([2*find(keep)-1,2*find(keep)].',[],1);(21:P.nLambda).'];
C.retainedZ=[reshape([2*find(keep)-1,2*find(keep)].',[],1);(21:25).'];
C.FstarResidual=norm(P.B*P.ystar-P.p,'fro');
C.targetSlack=min(-b(C.retainedRows));C.targetHC=aux.hC;C.targetLambda2=aux.lambda2;
C.targetMargin=[min(P.K1*aux.hInt(keep)-P.DeltaInt),min(P.K1*aux.hObs-P.DeltaObs),P.K1*aux.hC-P.DeltaC];
C.targetPatternGap=min(abs(dist-P.dc));C.targetUNorm=norm(un);
% Single-edge radial derivative bounds. Sum 125 spanning-tree products.
Ge=8/(exp(1)*P.dc);
He=(8/exp(1)+4*(256/exp(3)+54/exp(2)))/P.dc^2;
sm=(1+sqrt(3))/2;
C.Gh=625/512*(2*sqrt(5)/P.dc)*sqrt(sm*(1+sm)^3)*exp(-sm);
C.Hh=625/512*(8*He+24*Ge^2);
% Construct Theta0 as a Euclidean slow-state ball and shrink until the
% explicit row, graph and extended-safety inequalities hold uniformly.
r=2e-4;
for attempt=1:60
 dmax=max(dist(keep))+sqrt(2)*r;
 hmin=min(aux.hInt(keep));
 obsDist=max(vecnorm(reshape(P.ystar,[5,1,2])-reshape(P.obs,[1,P.nObs,2]),2,3),[],'all')+r;
 pairSlack=.5*(P.K1*(hmin-2*dmax*sqrt(2)*r)-P.DeltaInt)-2*dmax*(P.gamma+8*P.c+C.tildeL+P.K2)*r-2*r^2;
 chiMin=exp(-dmax^2/(P.dc^2-dmax^2));
 obsSlack=P.K1*(min(aux.hObs)-2*obsDist*r)-P.DeltaObs-2*obsDist*(P.gamma+8*P.c+C.tildeL+P.K2)*r-2*r^2;
 connSlack=(P.K1*(aux.hC-C.Gh*r)-P.DeltaC)/5-C.Gh*(P.gamma+8*P.c+C.tildeL+P.K2)*r-C.Hh*((1+8*P.c)*r)^2/5;
 C.slackLower=[P.barrierScale*chiMin*pairSlack,P.barrierScale*obsSlack,P.connScale*connSlack,P.deltaP];
 C.localBarrierLower=[hmin-2*dmax*sqrt(2)*r,min(aux.hObs)-2*obsDist*r,aux.hC-C.Gh*r];
 C.localExtendedLower=P.alpha1*C.localBarrierLower-[2*dmax*sqrt(2)*r,2*obsDist*r,C.Gh*r];
 C.localGraphLower=aux.lambda2-2*P.m*sqrt(2)*Ge*r;
 if all(C.slackLower>0) && min(C.localExtendedLower)>0 && C.localGraphLower>P.lambdaLower && sqrt(2)*r<C.targetPatternGap && r<P.nu/2,break;end
 r=r/2;
end
C.rTheta0=r;C.etaI=min(C.slackLower);
C.rowNormUpper=[sqrt((P.barrierScale*2*dmax)^2+2),P.barrierScale*2*obsDist,sqrt((P.connScale*C.Gh)^2+8^2),0];
C.Bbar=sqrt(sum([2*nnz(keep),5*P.nObs,5,5].*C.rowNormUpper.^2));
C.rLoc=min(1e-3,.5*min((C.slackLower(1:3)-C.etaI/2)./C.rowNormUpper(1:3)));
C.tubeSlackLower=C.slackLower-C.rowNormUpper*C.rLoc;
h0=min(1,P.xi);C.kappa=2+(C.Bbar^2/h0+1)*C.rLoc/C.etaI;
C.cf=min(h0/2,C.kappa*C.etaI/(2*C.rLoc)-C.Bbar^2/(2*h0));
C.mfUnder=.5*min(1,C.kappa);C.mfBar=.5*max(1,C.kappa);
% L_z bounds the field measured from equilibrium, since L(y)1=0.
C.Lz=sqrt((1+8*P.c)^2+(P.gamma+8*P.c+C.tildeL)^2);C.Lu=1;
C.sigma=(1+C.k)/(C.tildeL*C.Lz);
C.epsUpper=4*C.cf*C.alphaNom*C.sigma/((1+C.k+C.sigma*C.tildeL*C.Lz)^2+4*C.alphaNom*C.sigma*C.tildeL);
C.Q=[C.alphaNom,-.5*(1+C.k+C.sigma*C.tildeL*C.Lz);-.5*(1+C.k+C.sigma*C.tildeL*C.Lz),C.sigma*(C.cf/P.epsilon-C.tildeL)];
C.minEigQ=min(eig(C.Q));
C.rho=.5*min(C.mUnder*r^2,C.sigma*C.mfUnder*C.rLoc^2);
C.rTheta1=sqrt(C.rho/(2*C.mBar));C.rFast=sqrt(C.rho/(2*C.sigma*C.mfBar));
C.mCompositeUnder=min(C.mUnder,C.sigma*C.mfUnder);C.mCompositeBar=max(C.mBar,C.sigma*C.mfBar);
C.rateLower=C.minEigQ/(2*C.mCompositeBar);
C.analyticLocalPass=C.theorem1 && min(C.targetMargin)>0 && min(C.slackLower)>0 && min(C.localExtendedLower)>0 && C.localGraphLower>P.lambdaLower && min(C.tubeSlackLower)>=C.etaI/2 && P.epsilon<C.epsUpper && C.minEigQ>0 && C.rFast<C.rLoc;
C.regionalUniformCertificate=false;
C.scope='Explicit conservative local inequalities; floating-point evaluation; no interval certificate and no regional entrance claim.';
end
