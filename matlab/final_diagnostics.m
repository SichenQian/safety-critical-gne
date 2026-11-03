function final_diagnostics(name)
root=paper_data_root();folder=fullfile(root,'results',name);
D=load(fullfile(folder,'trajectory.mat'),'P','S');P=D.P;S=D.S;F=load(fullfile(folder,'final_full.mat'),'T','weights','dists');T=F.T;
active=F.dists<P.dc;ww=F.weights;wa=ww;wa(~active)=Inf;
tree=zeros(height(T),1);
for k=1:height(T)
 [sw,ord]=sort(ww(k,:),'descend');lab=1:5;count=0;
 for j=1:10
  e=ord(j);if ~active(k,e),continue;end
  a=lab(P.edges(e,1));b=lab(P.edges(e,2));if a==b,continue;end
  lab(lab==b)=a;count=count+1;
  if count==4,tree(k)=sw(j);break;end
 end
end
% Refine the minimum maximum-bottleneck spanning-tree weight on dense output.
[tm,treeMin]=minimize(@(t)metric(t,4),T.time,tree);
G=struct('maxComponents',max(T.components),'minimumActiveWeightComputed',min(wa,[],'all'),...
 'activeWeightBoundaryInfimum',0,'minimumBottleneckTreeWeight',treeMin,'minimumBottleneckTime',tm,...
 'underflowActiveSamples',nnz(T.underflowActiveCount),'stateCount',height(T));
if isempty(S.sensingCrossings),G.activeWeightBoundaryInfimum=NaN;end
% Barrier minima identify the exact geometry used by figure insets.
E=struct();fields={'collision','obstacle','tree'};columns={'hInt','hObs','hC'};
for j=1:3
 [t,h]=minimize(@(tt)metric(tt,j),T.time,T.(columns{j}));z=state(t);Y=positions(z);V=reshape(z(P.iV),2,5).';
 e=struct('time',t,'h',h,'Y',Y,'V',V);
 if j==1
  hs=sum((Y(P.edges(:,1),:)-Y(P.edges(:,2),:)).^2,2)-P.dsafe^2;[~,id]=min(hs);e.pair=P.edges(id,:);e.distance=norm(Y(e.pair(1),:)-Y(e.pair(2),:));e.geometricGap=e.distance-P.dsafe;
 elseif j==2
  hs=sum((reshape(Y,[5 1 2])-reshape(P.obs,[1 P.nObs 2])).^2,3)-reshape(P.rObs.^2,[1 P.nObs]);[~,id]=min(hs,[],'all','linear');[i,o]=ind2sub(size(hs),id);
  e.agent=i;e.obstacle=o;e.distance=norm(Y(i,:)-P.obs(o,:));e.geometricGap=e.distance-P.rObs(o);
 end
 E.(fields{j})=e;
end
% Bounds are a box INTERSECTED with a pairwise diameter and speed constraint.
cols=reshape([1:12:49;2:12:50],1,10);Y=[S.z(:,cols);S.internalZ(:,cols);S.denseZ(:,cols)];
region=struct('box',P.regionBox,'rectangleDiameter',norm(diff(P.regionBox,1,2)),...
 'pairwiseDiameterBound',P.regionDiameter,'collisionCoefficientBound',2*P.regionDiameter,...
 'observedPairwiseDiameter',max(T.maxDistance),'observedIndividualSpeed',max(T.maxSpeed),...
 'xRange',[min(Y(:,1:2:9),[],'all') max(Y(:,1:2:9),[],'all')],'yRange',[min(Y(:,2:2:10),[],'all') max(Y(:,2:2:10),[],'all')]);
region.trajectoryWithinRegion=region.observedPairwiseDiameter<=P.regionDiameter && region.observedIndividualSpeed<=P.regionSpeedMax && region.xRange(1)>=P.regionBox(1,1)&&region.xRange(2)<=P.regionBox(1,2)&&region.yRange(1)>=P.regionBox(2,1)&&region.yRange(2)<=P.regionBox(2,2);
write_json(fullfile(folder,'graph_statistics.json'),G);write_json(fullfile(folder,'geometry_extrema.json'),E);write_json(fullfile(folder,'region_check.json'),region);
save(fullfile(folder,'figure_diagnostics.mat'),'G','E','region','tree','-v7');
fprintf('FINAL_GRAPH_GEOMETRY_DONE %s\n',name);
 function z=state(t)
  ends=cellfun(@(s)s.x(end),S.solutions);id=find(ends>=t,1);z=deval(S.solutions{id},t);
 end
 function Y=positions(z)
  Y=zeros(5,2);for posAgent=1:5,posIndex=(posAgent-1)*10+2*posAgent-1;Y(posAgent,:)=z(posIndex:posIndex+1).';end
 end
 function v=metric(t,j)
  z=state(t);Y=positions(z);
  if j==1,v=min(sum((Y(P.edges(:,1),:)-Y(P.edges(:,2),:)).^2,2)-P.dsafe^2);
  elseif j==2,v=min(sum((reshape(Y,[5 1 2])-reshape(P.obs,[1 P.nObs 2])).^2,3)-reshape(P.rObs.^2,[1 P.nObs]),[],'all');
  elseif j==3,v=audit_core('tree',P,Y,zeros(5,2));
  else
   [W,~]=audit_core('graph',P,Y);weights=W(sub2ind([5 5],P.edges(:,1),P.edges(:,2)));[sw,ord]=sort(weights,'descend');lab=1:5;count=0;v=0;
   for jj=1:10,e=ord(jj);a=lab(P.edges(e,1));b=lab(P.edges(e,2));if a==b,continue;end;lab(lab==b)=a;count=count+1;if count==4,v=sw(jj);break;end,end
  end
 end
end
function [tm,v]=minimize(fun,t,y)
[v,k]=min(y);tm=t(k);idx=find(y(2:end-1)<=y(1:end-2)&y(2:end-1)<=y(3:end))+1;[~,ord]=sort(y(idx));idx=idx(ord(1:min(16,numel(ord))));
for k=idx(:).', [tt,f]=fminbnd(fun,t(k-1),t(k+1),optimset('TolX',1e-10,'Display','off'));if f<v,v=f;tm=tt;end,end
end
