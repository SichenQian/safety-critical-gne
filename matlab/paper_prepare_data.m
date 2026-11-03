function D=paper_prepare_data(source)
% Build every plotted quantity from the six saved numerical experiments.
if nargin<1,source=paper_data_root();end
A=load(fullfile(source,'results','selected','trajectory.mat'),'P','S');
B=load(fullfile(source,'results','untightened','trajectory.mat'),'S');
P=A.P;E=jsondecode(fileread(fullfile(source,'results','untightened','geometry_extrema.json')));
assert(isequal(A.S.z(1,:),B.S.z(1,:)));
cols=reshape([1:12:49;2:12:50],1,10);
D=struct('P',P,'E',E,'timeT',A.S.t,'timeU',B.S.t,'YT',A.S.z(:,cols),'YU',B.S.z(:,cols));
F=load(fullfile(source,'results','selected','final_full.mat'),'T');
G=load(fullfile(source,'results','untightened','final_full.mat'),'T');
D.T=thin(F.T);D.U=thin(G.T);
D.zoomSafety=G.T(G.T.time>=55&G.T.time<=140,{'time','hInt','hObs'});
D.zoomTree=G.T(abs(G.T.time-E.tree.time)<=40,{'time','hC'});
R=jsondecode(fileread(fullfile(source,'results','selected','final_full.json')));D.events=R.crossings;
D.obstacle=geometry(E.obstacle,'obstacle');D.collision=geometry(E.collision,'collision');
L=load(fullfile(source,'results','local','final_full.mat'),'T');
M=load(fullfile(source,'results','local','trajectory.mat'),'P','C','initial');
mask=L.T.time<=5*P.epsilon/P.xi;
D.local=struct('time',L.T.time,'slow',L.T.slowError,'fast',L.T.fastError,...
 'fastTime',L.T.time(mask)/P.epsilon,'fastLayer',L.T.fastError(mask),...
 'initialState',M.P.initialState,'initial',M.initial,'certificate',M.C);
D.slowFit=polyfit(L.T.time(L.T.time>=20&L.T.time<=120),log(L.T.slowError(L.T.time>=20&L.T.time<=120)),1);
D.fastFit=polyfit(L.T.time(mask),log(L.T.fastError(mask)),1);
D.description='Data-derived trajectories, physical margins, exact support events and independent local error norms.';
 function Z=geometry(e,kind)
  if strcmp(kind,'obstacle'),q=e.Y(e.agent,:)-P.obs(e.obstacle,:);v=e.V(e.agent,:);rad=P.rObs(e.obstacle);
  else,q=e.Y(e.pair(2),:)-e.Y(e.pair(1),:);v=e.V(e.pair(2),:)-e.V(e.pair(1),:);rad=P.dsafe;end
  width=1.8*abs(e.geometricGap);ex=floor(log10(width));scale=10^ex;
  normal=q/norm(q);tan=[-normal(2) normal(1)];
  dt=min(.5,2*width/max(abs(v*tan.'),1e-8));tt=linspace(e.time-dt,e.time+dt,801);zz=states(B.S,tt);
  if strcmp(kind,'obstacle'),ii=(e.agent-1)*10+2*e.agent-1;qq=zz(ii:ii+1,:).'-P.obs(e.obstacle,:);
  else,ii=(e.pair(1)-1)*10+2*e.pair(1)-1;jj=(e.pair(2)-1)*10+2*e.pair(2)-1;qq=(zz(jj:jj+1,:)-zz(ii:ii+1,:)).';end
  xx=linspace(-width,width,501);yy=-xx.^2./(sqrt(rad^2-xx.^2)+rad);
  Z=struct('time',tt,'trace',(qq-rad*normal)*[tan.' normal.']/scale,...
   'boundary',[xx;yy].'/scale,'width',width/scale,'exponent',ex,'gap',e.geometricGap/scale);
 end
end
function T=thin(T)
ix=unique([1:15:height(T),height(T)]);
for c={'hInt','hObs','hC','lambda2','tracking','deltaBall'}
 [~,a]=min(T.(c{1}));[~,b]=max(T.(c{1}));ix=union(ix,[a b]);
end
T=T(ix,:);
end
function z=states(S,t)
ends=cellfun(@(s)s.x(end),S.solutions);z=zeros(size(S.z,2),numel(t));
for k=1:numel(t),j=find(ends>=t(k),1);z(:,k)=deval(S.solutions{j},t(k));end
end
