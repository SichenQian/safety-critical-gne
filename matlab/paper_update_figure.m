function paper_update_figure(f,name,D)
% Populate final layout templates exclusively from the numerical data.
switch name
 case 'trajectory'
  for i=1:5
   line(f,sprintf('trajectory_T_%d',i),D.YT(:,2*i-1),D.YT(:,2*i));
   line(f,sprintf('trajectory_U_%d',i),D.YU(:,2*i-1),D.YU(:,2*i));
   line(f,sprintf('paper_initial_%d',i),D.P.y0(i,1),D.P.y0(i,2));
   line(f,sprintf('paper_target_%d',i),D.P.ystar(i,1),D.P.ystar(i,2));
  end
  meta=f.UserData.paperBinding;
  for c={'obstacle','collision'}
   kind=c{1};Z=D.(kind);assert(Z.exponent==meta.([kind 'Exponent']));
   line(f,['paper_' kind '_trace'],Z.trace(:,1),Z.trace(:,2));
   line(f,['paper_' kind '_boundary'],Z.boundary(:,1),Z.boundary(:,2));
   line(f,['paper_' kind '_gap'],0,Z.gap);
   h=findall(f,'Tag',['paper_' kind '_forbidden']);x=Z.boundary(:,1).';y=Z.boundary(:,2).';
   set(h,'XData',[x fliplr(x)],'YData',[y -Z.width*ones(size(x))]);
   if strcmp(kind,'obstacle'),a=findall(f,'Tag','zoom_A');else,a=findall(f,'Tag','zoom_B');end
   a.XLim=[-Z.width Z.width];a.YLim=[-Z.width Z.width];
  end
  q=D.E.obstacle.Y(D.E.obstacle.agent,:);pp=D.E.collision.Y(D.E.collision.pair,:);
  line(f,'paper_obstacle_point',q(1),q(2));line(f,'paper_collision_point_1',pp(1,1),pp(1,2));line(f,'paper_collision_point_2',pp(2,1),pp(2,2));
  ra=findall(f,'Tag','region_A');ra.Position=[q+meta.regionAOffset meta.regionASize];
  rb=findall(f,'Tag','region_B');lo=min(pp,[],1)-meta.regionBLowerPad;hi=max(pp,[],1)+meta.regionBUpperPad;rb.Position=[lo hi-lo];
  for k={'A','B'}
   key=k{1};h=findall(f,'Tag',['paper_leader_' key]);if isempty(h),continue;end
   if key=='A',box=ra.Position;else,box=rb.Position;end
   idx=meta.(['leader' key 'Index']);tip=box(1:2)+meta.(['leader' key 'Fraction']).*box(3:4);
   x=h.XData;y=h.YData;x(idx)=tip(1);y(idx)=tip(2);set(h,'XData',x,'YData',y);
  end
  for i=1:D.P.nObs
   h=findall(f,'Tag',sprintf('paper_static_obstacle_%d',i));
   th=linspace(0,2*pi,401);set(h,'XData',D.P.obs(i,1)+D.P.rObs(i)*cos(th),'YData',D.P.obs(i,2)+D.P.rObs(i)*sin(th));
  end
 case 'safety'
  a=findall(f,'Tag','main_a');z=findall(f,'Tag','zoom_a');
  line(a,'collision_tightened',D.T.time,D.T.hInt);line(a,'collision_untightened',D.U.time,D.U.hInt);
  line(a,'obstacle_tightened',D.T.time,D.T.hObs);line(a,'obstacle_untightened',D.U.time,D.U.hObs);
  line(z,'collision_untightened',D.zoomSafety.time,D.zoomSafety.hInt/1e-4);
  line(z,'obstacle_untightened',D.zoomSafety.time,D.zoomSafety.hObs/1e-4);
 case 'connectivity'
  a=findall(f,'Tag','main_b');z=findall(f,'Tag','zoom_b');
  line(a,'certificate_tightened',D.T.time,D.T.hC);line(a,'certificate_untightened',D.U.time,D.U.hC);
  line(a,'lambda2_tightened',D.T.time,D.T.lambda2);line(a,'lambda2_untightened',D.U.time,D.U.lambda2);
  line(z,'certificate_zoom_untightened',D.zoomTree.time,D.zoomTree.hC/1e-3);
 case 'edge'
  assert(size(D.events,1)==2&&isequal(D.events(:,2:3),[2 5;3 4]));
  a=findall(f,'Tag','main_c');t=[0 D.events(1,1) D.events(1,1) D.events(2,1) D.events(2,1) a.XLim(2)];
  line(f,'edge_2_5_active',t,[0 0 1 1 1 1]);line(f,'edge_3_4_active',t,[1 1 1 1 0 0]);
  line(f,'total_active_edges',t,[9 9 10 10 9 9]);
 case 'exp'
  line(f,'slow_error',D.local.time,D.local.slow);line(f,'fast_error',D.local.time,D.local.fast);
  line(f,'fast_error_inset',D.local.fastTime,D.local.fastLayer);
 otherwise,error('Unknown figure %s',name);
end
end
function line(parent,tag,x,y)
h=findall(parent,'Type','line','Tag',tag);
if isempty(h)&&(contains(tag,'_point')||endsWith(tag,'_gap')),return;end
assert(numel(h)==1,'Missing/ambiguous data binding: %s',tag);
set(h,'XData',x,'YData',y);
end
