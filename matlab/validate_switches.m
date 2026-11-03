function validate_switches(name)
if nargin<1,name='selected';end
root=paper_data_root();folder=fullfile(root,'results',name);
D=load(fullfile(folder,'trajectory.mat'));P=D.P;S=D.S;rows=[];
for k=1:size(S.sensingCrossings,1)
 ev=S.sensingCrossings(k,:);tt=ev(1)+[-.01 -.001 -1e-5 0 1e-5 .001 .01];
 for t=tt
  ends=cellfun(@(s)s.x(end),S.solutions);j=find(ends>=t,1);z=deval(S.solutions{j},t);A=audit_state(P,z,true);
  [M,b,Mc,bc]=normalized_reference_rows(P,z,true,A);
  Q=audit_qp(diag([ones(10,1);P.xi*ones(25,1)]),[-A.nom;zeros(25,1)],M,b);
  C=audit_qp(eye(10),-A.nom,Mc,bc);
  e=find(all(P.edges==ev(2:3),2));gap=norm(A.Y(ev(2),:)-A.Y(ev(3),:))-P.dc;
  rows(end+1,:)=[k,t,ev(2:3),gap,gap<0,A.hc,A.lambda2,Q.referenceValid,C.referenceValid,norm(z(P.iU)-Q.w(1:10)),A.deltaBall,A.trueQmax,A.psiPair(e)];
 end
end
T=array2table(rows,'VariableNames',{'event','time','i','j','distanceMinusDc','active','hC','lambda2','QPvalid','centralQPvalid','tracking','deltaBall','qInt','qObs','qC','psiPair'});
writetable(T,fullfile(folder,'switch_neighborhoods.csv'));
assert(all(T.QPvalid)&all(T.centralQPvalid)&all(T.hC>0)&all(T.lambda2>P.lambdaLower)&all(T.tracking<P.trackingBudget)&all(T.deltaBall<P.connMismatchBudget));
assert(all(T.qInt<=0)&all(T.qObs<=0)&all(T.qC<=0)&all(T.psiPair>0));
for k=1:size(S.sensingCrossings,1)
 ix=find(T.event==k);assert(T.distanceMinusDc(ix(1))*T.distanceMinusDc(ix(end))<0);
end
write_json(fullfile(folder,'switch_validation.json'),struct('events',S.sensingCrossings,'neighborhoodChecks',height(T),'allPassed',true,'supportDefinition','distance < dc; not floating-point positive weight','continuousTimeProof',false));
end
