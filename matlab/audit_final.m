function R=audit_final(name,mode)
% Dense, internal-step and local-extremum checks. All are numerical checks.
% A physical or tracking supremum over an entire region is never inferred.
if nargin<2,mode='full';end
root=paper_data_root();folder=fullfile(root,'results',name);
tmp=load(fullfile(folder,'trajectory.mat'),'P','S');P=tmp.P;S=tmp.S;
if strcmp(mode,'quick'),t=S.t;Z=S.z;elseif strcmp(mode,'preflight')
 [t,ix]=unique([S.t;S.internalT]);ZZ=[S.z;S.internalZ];Z=ZZ(ix,:);
else
 [t,ix]=unique([S.t;S.internalT;S.denseT]);ZZ=[S.z;S.internalZ;S.denseZ];Z=ZZ(ix,:);
end
nr=numel(t);data=zeros(nr,30);dists=zeros(nr,10);weights=dists;
H=diag([ones(10,1);P.xi*ones(25,1)]);
for k=1:nr
 if mod(k-1,10000)==0,fprintf('FINAL %s %d/%d\n',name,k,nr);end
 z=Z(k,:).';D=audit_state(P,z,S.tight);dist=vecnorm(D.Y(P.edges(:,1),:)-D.Y(P.edges(:,2),:),2,2);active=dist<P.dc;
 [A,b,Ac,bc,meta]=normalized_reference_rows(P,z,S.tight,D);
 Q=audit_qp(H,[-D.nom;zeros(25,1)],A,b);
 % Evaluate the decomposed input as a candidate for the PHYSICAL centralized
 % QP. Estimated connectivity rows need not sum to the physical row: a
 % positive value here does not prove centralized-QP infeasibility.
 track=NaN;if Q.referenceValid,track=norm(z(P.iU)-Q.w(1:10));end
 q=D.q;qp=q(1:10);qp(~active)=-Inf;
 pbar=dist.^2-P.dsafe^2;
 obar=zeros(5*P.nObs,1);for o=1:P.nObs,obar((o-1)*5+(1:5))=sum((D.Y-P.obs(o,:)).^2,2)-P.rObs(o)^2;end
 centralWitness=NaN;if Q.referenceValid,centralWitness=max(Ac*Q.w(1:10)+bc);end
 data(k,:)=[t(k),min(pbar),min(obar),D.hc,D.lambda2,max(qp),max(q(P.cObs)),q(P.cConn),track,D.deltaActual,D.deltaBall,...
  D.estimateError,D.slowError,D.fastLocalError,max(0,Q.primal),Q.stationarity,Q.complementarity,Q.dualViolation,...
  Q.referenceValid,max(vecnorm(D.V,2,2)),max(dist),norm(z(P.iU)),norm(Q.w(1:10)),nnz(active),meta.underflowActiveCount,...
  min([D.psiPair(active);D.psiObs;D.psiC]),D.components,centralWitness,min(z(P.iL)),max(abs(D.Y),[],'all')];
 dists(k,:)=dist;weights(k,:)=D.weights;
end
names={'time','hInt','hObs','hC','lambda2','qInt','qObs','qC','tracking','deltaActual','deltaBall','estimateError','slowError','fastError','qpPrimal','qpStationarity','qpComplementarity','qpDual','qpValid','maxSpeed','maxDistance','inputNorm','optimizerInputNorm','edgeCount','underflowActiveCount','minPsi','components','centralWitnessResidual','minDual','maxPositionAbs'};
T=array2table(data,'VariableNames',names);writetable(T,fullfile(folder,['final_' mode '.csv']));
writetable(array2table([t,dists,weights],'VariableNames',[{'time'},compose('distance_%d',1:10),compose('weight_%d',1:10)]),fullfile(folder,['topology_' mode '.csv']));
R=struct('samples',nr,'validQP',nnz(T.qpValid),'minimumBarriers',min(data(:,2:4),[],1),'minimumLambda2',min(T.lambda2),...
 'maxECBF',max(data(:,6:8),[],1),'maxTracking',max(T.tracking),'maxDeltaActual',max(T.deltaActual),'maxDeltaBall',max(T.deltaBall),...
 'maxKKT',max(data(:,15:18),[],1),'maxSpeed',max(T.maxSpeed),'maxDistance',max(T.maxDistance),'maxInputNorm',max(T.inputNorm),...
 'maxOptimizerInputNorm',max(T.optimizerInputNorm),'edgeCountRange',[min(T.edgeCount),max(T.edgeCount)],...
 'underflowSamples',nnz(T.underflowActiveCount),'minimumPsi',min(T.minPsi),'maxComponents',max(T.components),...
 'maxCentralWitnessResidual',max(T.centralWitnessResidual),'minimumDual',min(T.minDual),'continuousTimeProof',false,...
 'KKTcoordinates','Positive collision weights analytically cancelled for distance<dc; zero rows retained at distance>=dc.',...
 'crossings',S.sensingCrossings);
write_json(fullfile(folder,['final_' mode '.json']),R);
save(fullfile(folder,['final_' mode '.mat']),'P','R','T','dists','weights','-v7');
% Find local maxima by bounded scalar optimization of the saved interpolant.
if ~strcmp(mode,'quick') && ~isempty(S.solutions)
 peaks=struct();fields={'tracking','deltaBall'};
 for j=1:2
  val=T.(fields{j});idx=find(val(2:end-1)>=val(1:end-2) & val(2:end-1)>=val(3:end))+1;
  [~,ord]=sort(val(idx),'descend');idx=idx(ord(1:min(12,numel(ord))));
  out=zeros(0,3);
  for k=idx(:).'
   lo=t(k-1);hi=t(k+1);
   [tp,fv]=fminbnd(@(tt)-metric(tt,j),lo,hi,optimset('TolX',1e-9,'Display','off'));
   out(end+1,:)=[tp,-fv,val(k)]; %#ok<AGROW>
  end
  peaks.(fields{j})=out;
 end
 peaks.scope='Local maximization on numerical ode15s dense interpolants; no global or interval-arithmetic proof.';
 if strcmp(mode,'full'),peakFile='refined_peaks.json';else,peakFile=['refined_peaks_' mode '.json'];end
 write_json(fullfile(folder,peakFile),peaks);
end
    function v=metric(tt,j)
      ends=cellfun(@(s)s.x(end),S.solutions);id=find(ends>=tt,1);zz=deval(S.solutions{id},tt);
      if isfield(S,'solutionCoordinates') && strcmp(S.solutionCoordinates,'equilibrium_error'),zz=zz+S.equilibrium;end
      D=audit_state(P,zz,S.tight);
      if j==2,v=D.deltaBall;else
        [A,b]=normalized_reference_rows(P,zz,S.tight,D);Q=audit_qp(H,[-D.nom;zeros(25,1)],A,b);
        if ~Q.referenceValid,error('Unresolved QP in dense peak optimization');end
        v=norm(zz(P.iU)-Q.w(1:10));
      end
    end
end
