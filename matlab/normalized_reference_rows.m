function [A,b,Ac,bc,meta]=normalized_reference_rows(P,z,tight,D)
% OFFLINE ONLY. Cancel strictly positive chi analytically before evaluating.
% This prevents exponential underflow from deleting a mathematically active
% collision row. Exact inactivity is determined by distance >= dc.
% KKT multipliers refer to these row coordinates, not the raw chi-scaled row.
A=D.A;b=D.b;Ac=D.Ac;bc=D.bc;
X=reshape(z(P.iX),10,5).';[~,~,~,aux]=audit_core('rows',P,X,D.V,tight);
dist=vecnorm(D.Y(P.edges(:,1),:)-D.Y(P.edges(:,2),:),2,2);active=dist<P.dc;
for e=1:P.m
 i=P.edges(e,1);j=P.edges(e,2);d=D.Y(i,:)-D.Y(j,:);v=D.V(i,:)-D.V(j,:);h=sum(d.^2)-P.dsafe^2;
 rr=2*e-1:2*e;A(rr,:)=0;b(rr)=0;Ac(e,:)=0;bc(e)=0;
 if ~active(e),continue;end
 for side=1:2
  if side==1,a=i;sgn=1;own=2*e-1;other=2*e;else,a=j;sgn=-1;own=2*e;other=2*e-1;end
  g=sgn*2*d;r=rr(side);A(r,2*a-1:2*a)=-P.barrierScale*g;A(r,10+own)=1;A(r,10+other)=-1;
  b(r)=P.barrierScale*(-g*aux.f0(a,:).'-P.K2*g*D.V(a,:).'-sum(v.^2)-.5*P.K1*h+.5*tight*parameter_at(P.DeltaInt,e));
 end
 Ac(e,2*i-1:2*i)=-2*d;Ac(e,2*j-1:2*j)=2*d;
 bc(e)=-2*d*(aux.f0(i,:)-aux.f0(j,:)).'-2*sum(v.^2)-2*P.K2*d*v.'-P.K1*h+tight*parameter_at(P.DeltaInt,e);
end
meta=struct('active',active,'distances',dist,'underflowActiveCount',nnz(active & D.weights==0));
end
