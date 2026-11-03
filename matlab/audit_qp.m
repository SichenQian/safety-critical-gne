function Q=audit_qp(H,f,A,b)
% Offline only. Positive row normalization, multiplier mapped back to raw rows.
norms=vecnorm(A,2,2);keep=norms>0;scale=ones(size(b));scale(keep)=min(1./norms(keep),100);
% Strictly positive sensing weights must not make collision rows disappear
% inside an absolute QP tolerance. Remove their positive row factor even
% when chi is tiny. Keep the cap for vanishing dissipativity rows.
pairCount=0;if numel(f)==35 && numel(b)>=35,pairCount=20;elseif numel(f)==10 && numel(b)>=21,pairCount=10;end
pairRows=(1:pairCount).';pairRows=pairRows(keep(pairRows));scale(pairRows)=min(1./norms(pairRows),realmax);
As=A.*scale;bs=b.*scale;opt=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex','ConstraintTolerance',1e-11,'OptimalityTolerance',1e-11,'MaxIterations',2000);
T=ones(size(f));if numel(f)==35,T(11:30)=sqrt(.06);T(31:35)=sqrt(.1);end
wt=-H\f;method='interior_point_bounded_scaling';
if max(A*wt+b)<=0
 w=wt;flag=1;output.iterations=0;lm.ineqlin=zeros(size(b));method='exact_unconstrained_feasible';
else
 [ww,~,flag,output,lm]=quadprog((H.*T).*T.',f.*T,As.*T.',-bs,[],[],[],[],[],opt);w=T.*ww;
end
rawflag=flag;rawResidual=NaN;if numel(w)==numel(f),rawResidual=max(A*w+b);end
% Active-set KKT polishing, accepted only after checking every original row.
if numel(w)==numel(f) && all(isfinite(w)) && numel(lm.ineqlin)==numel(b)
 ia=find(lm.ineqlin>1e-8 & abs(As*w+bs)<1e-5);
 if ~isempty(ia)
  K=[H,As(ia,:).';As(ia,:),zeros(numel(ia))];rr=[-f;-bs(ia)];
  if rcond(K)>1e-13
   sol=K\rr;wp=sol(1:numel(f));lp=sol(numel(f)+1:end);
   if min(lp)>=-1e-10 && max(As*wp+bs)<1e-9 && norm(H*wp+f+As(ia,:).'*lp,inf)<1e-9
    w=wp;lm.ineqlin=zeros(size(b));lm.ineqlin(ia)=max(lp,0);flag=1;method='verified_active_KKT_polish';
   end
  end
 end
end
% If one conditioning choice fails, retry the SAME QP in the original
% package's fixed row coordinates. Feasible/KKT checks, not exitflags alone,
% decide whether the reference can be used. No online control is changed.
good=false;
if numel(w)==numel(f) && all(isfinite(w)) && numel(lm.ineqlin)==numel(b)
 lr=scale.*lm.ineqlin;rr=A*w+b;
 good=flag>0 && max(As*w+bs)<=1e-8 && norm(H*w+f+A.'*lr,inf)<1e-8 && norm(lr.*rr,inf)<1e-8;
end
if ~good
 sr=ones(size(b));if numel(f)==35 && numel(b)>=35,sr=[ones(numel(b)-10,1)/sqrt(.06);ones(5,1)/sqrt(.1);ones(5,1)];end
 [wr,~,fr,~,lr]=quadprog((H.*T).*T.',f.*T,(A.*sr).*T.',-b.*sr,[],[],[],[],wt./T,opt);
 if numel(wr)==numel(f) && all(isfinite(wr)) && numel(lr.ineqlin)==numel(b)
  wr=T.*wr;lamr=sr.*lr.ineqlin;rr=A*wr+b;
  if fr>0 && max(As*wr+bs)<=1e-8 && norm(H*wr+f+A.'*lamr,inf)<1e-7 && norm(lamr.*rr,inf)<1e-7
   w=wr;lm.ineqlin=lamr./scale;flag=fr;method='verified_fixed_scaling_retry';
  end
 end
end
Q=struct('exitflag',flag,'feasible',false,'w',nan(size(f)),'lambda',nan(size(b)),...
 'primal',NaN,'stationarity',NaN,'complementarity',NaN,'dualViolation',NaN,'normalizedPrimal',NaN,...
 'normalizedMinSlack',NaN,'minPositiveRowNorm',NaN,'conditionActive',NaN,'iterations',output.iterations,'conflict',struct(),...
 'rawSolverExitflag',rawflag,'rawSolverPrimal',rawResidual,'method',method,'referenceValid',false,'classification','solver_unresolved',...
 'feasibleWitnessPrimal',NaN,'feasibleWitnessNormalizedPrimal',NaN);
if any(keep),Q.minPositiveRowNorm=min(norms(keep));end
if numel(w)==numel(f) && all(isfinite(w)) && numel(lm.ineqlin)==numel(b)
 lam=scale.*lm.ineqlin;r=A*w+b;rn=As*w+bs;Q.w=w;Q.lambda=lam;
 Q.primal=max(r);Q.normalizedPrimal=max(rn);Q.stationarity=norm(H*w+f+A.'*lam,inf);
 Q.complementarity=norm(lam.*r,inf);Q.dualViolation=max(0,-min(lam));
 Q.feasible=Q.primal<=1e-8 && Q.normalizedPrimal<=1e-7;
 if Q.feasible,Q.feasibleWitnessPrimal=Q.primal;Q.feasibleWitnessNormalizedPrimal=Q.normalizedPrimal;end
 Q.referenceValid=Q.feasible && Q.stationarity<1e-7 && Q.complementarity<1e-7 && Q.dualViolation<1e-8;
 if Q.referenceValid,Q.classification='feasible_KKT_reference';elseif Q.feasible,Q.classification='feasible_point_optimizer_unresolved';end
 if any(keep),Q.normalizedMinSlack=min(-rn(keep));end
 ia=keep & rn>-1e-7;if any(ia),sv=svd(As(ia,:),'econ');Q.conditionActive=max(sv)/max(min(sv),realmin);end
end
if ~Q.feasible
 % Minimum normalized common relaxation; a failed numerical LP is not a proof.
 opts=optimoptions('linprog','Display','off','OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);
 [v,val,ef,~,dl]=linprog([zeros(size(f));1],[As,-ones(size(b))],-bs,[],[],[-inf(size(f));0],[],opts);
 Q.conflict=struct('relaxationExitflag',ef,'minimumNormalizedRelaxation',val,'rowIndices',[],'rawResidual',[],'normalizedResidual',[],'dualWeights',[]);
 if ef>0
  rr=As*v(1:end-1)+bs;ix=find(rr>max(1e-8,.9*val) | dl.ineqlin>1e-7);
  Q.conflict.rowIndices=ix;Q.conflict.rawResidual=A(ix,:)*v(1:end-1)+b(ix);Q.conflict.normalizedResidual=rr(ix);Q.conflict.dualWeights=dl.ineqlin(ix);
  Q.conflict.relaxedPoint=v(1:end-1);Q.conflict.dualBalance=norm(As.'*dl.ineqlin,inf);
  Q.conflict.farkasValue=bs.'*dl.ineqlin;Q.conflict.dualMinimum=min(dl.ineqlin);
  Q.conflict.selectedA=A(ix,:);Q.conflict.selectedB=b(ix);Q.conflict.rowNormalization=scale(ix);
  if val<=1e-8 && max(As*v(1:end-1)+bs)<=1e-7 && max(A*v(1:end-1)+b)<=1e-8
   Q.feasible=true;Q.classification='feasible_LP_witness_optimizer_unresolved';
   Q.feasibleWitnessPrimal=max(A*v(1:end-1)+b);Q.feasibleWitnessNormalizedPrimal=max(As*v(1:end-1)+bs);
   % A feasible LP point initializes a distinct active-set QP retry.
   oa=optimoptions('quadprog','Display','off','Algorithm','active-set','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',5000);
   [wa,~,fa,~,la]=quadprog((H.*T).*T.',f.*T,As.*T.',-bs,[],[],[],[],v(1:end-1)./T,oa);
   if numel(wa)==numel(f) && all(isfinite(wa)) && numel(la.ineqlin)==numel(b)
    wa=T.*wa;lr=scale.*la.ineqlin;rr=A*wa+b;st=norm(H*wa+f+A.'*lr,inf);cp=norm(lr.*rr,inf);
    Q.activeSetExitflag=fa;Q.activeSetKKT=[max(rr),st,cp];
    if max(rr)<=1e-8 && max(As*wa+bs)<1e-7 && st<1e-7 && cp<1e-7 && min(lr)>-1e-8
     Q.w=wa;Q.lambda=lr;Q.exitflag=fa;Q.primal=max(rr);Q.stationarity=st;Q.complementarity=cp;Q.dualViolation=max(0,-min(lr));Q.normalizedPrimal=max(As*wa+bs);
     rn=As*wa+bs;if any(keep),Q.normalizedMinSlack=min(-rn(keep));end
     ia=keep & rn>-1e-7;Q.conditionActive=NaN;if any(ia),sv=svd(As(ia,:),'econ');Q.conditionActive=max(sv)/max(min(sv),realmin);end
     Q.referenceValid=true;Q.classification='feasible_active_set_KKT_reference';Q.method='verified_LP_active_set_retry';
    end
   end
  elseif val>1e-8
   Q.classification='numerically_infeasible_positive_phaseI';
  end
 end
end
end
