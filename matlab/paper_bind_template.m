function paper_bind_template(f,name,D)
% Tag numerical objects in the final manually edited paper layout.
if ~strcmp(name,'trajectory'),return;end
a=findall(f,'Type','axes','Tag','trajectory');
meta=struct();
for pair={'obstacle','collision'}
 kind=pair{1};if strcmp(kind,'obstacle'),atag='zoom_A';else,atag='zoom_B';end
 z=findall(f,'Type','axes','Tag',atag);assert(numel(z)==1);
 for h=findall(z,'Type','line').'
  n=numel(h.XData);
  if n==801,suffix='trace';elseif n==501,suffix='boundary';elseif n==1,suffix='gap';else,error('Unexpected geometry line.');end
  h.Tag=['paper_' kind '_' suffix];
 end
 h=findall(z,'Type','patch');assert(numel(h)==1);h.Tag=['paper_' kind '_forbidden'];
 meta.([kind 'Exponent'])=D.(kind).exponent;
end
q=D.E.obstacle.Y(D.E.obstacle.agent,:);
pp=D.E.collision.Y(D.E.collision.pair,:);
for h=findall(a,'Type','line').'
 if ~isempty(h.Tag),continue;end
 if numel(h.XData)==1
  point=[h.XData h.YData];
  if norm(point-q)<1e-6,h.Tag='paper_obstacle_point';
  elseif norm(point-pp(1,:))<1e-6,h.Tag='paper_collision_point_1';
  elseif norm(point-pp(2,:))<1e-6,h.Tag='paper_collision_point_2';
  else
   [v,i]=min(vecnorm(D.P.y0-point,2,2));[w,j]=min(vecnorm(D.P.ystar-point,2,2));
   if v<1e-8,h.Tag=sprintf('paper_initial_%d',i);elseif w<1e-8,h.Tag=sprintf('paper_target_%d',j);end
  end
 end
end
ra=findall(a,'Tag','region_A');rb=findall(a,'Tag','region_B');assert(numel(ra)==1&&numel(rb)==1);
meta.regionAOffset=ra.Position(1:2)-q;meta.regionASize=ra.Position(3:4);
meta.regionBLowerPad=min(pp,[],1)-rb.Position(1:2);
meta.regionBUpperPad=rb.Position(1:2)+rb.Position(3:4)-max(pp,[],1);
for h=findall(a,'Type','line').'
 if numel(h.XData)~=2||~all(h.Color==0),continue;end
 points=[h.XData(:),h.YData(:)];
 da=vecnorm(points-(ra.Position(1:2)+ra.Position(3:4)/2),2,2);
 db=vecnorm(points-(rb.Position(1:2)+rb.Position(3:4)/2),2,2);
 if min(da)<min(db),key='A';box=ra.Position;[~,idx]=min(da);
 else,key='B';box=rb.Position;[~,idx]=min(db);end
 h.Tag=['paper_leader_' key];meta.(['leader' key 'Index'])=idx;
 meta.(['leader' key 'Fraction'])=(points(idx,:)-box(1:2))./box(3:4);
end
for h=findall(a,'Type','patch').'
 if numel(h.XData)~=401,continue;end
 point=[mean(h.XData) mean(h.YData)];[~,i]=min(vecnorm(D.P.obs-point,2,2));
 h.Tag=sprintf('paper_static_obstacle_%d',i);
end
f.UserData=struct('paperBinding',meta);
end
