function S=simulate_local_error(P,zstar)
% Exact translation of the SAME full projected ODE, including omitted rows.
% L(y)*(1 kron ystar)=0 and B*ystar=p remove cancellation near equilibrium.
% No quasi-steady substitution, no removal of rows and no online QP solves.
Pe=P;Pe.errorCoordinates=true;H=diag([ones(10,1);P.xi*ones(25,1)]);
q0=P.initialState-zstar;start=tic;
opts=odeset('RelTol',P.RelTol,'AbsTol',P.AbsTol,'MaxStep',P.maxStep,'Jacobian',@jac);
sol=ode15s(@rhs,[0 P.T],q0,opts);
tt=P.outputTimes(:);dt=unique([(0:P.denseDt:P.T).';tt]);
S=struct('t',tt,'z',deval(sol,tt).'+zstar.','tight',true,'wallSeconds',toc(start),'initialQPExit',1,'completed',true,...
 'faceEvents',zeros(0,2),'minDualNumerical',min(sol.y(P.iL,:),[],'all'),'internalT',sol.x(:),'internalZ',sol.y.'+zstar.',...
 'denseT',dt,'denseZ',deval(sol,dt).'+zstar.','sensingCrossings',zeros(0,7),'solver','ode15s in equilibrium error coordinates',...
 'controller','Original full projected primal-dual ODE; exact coordinate translation','solutionCoordinates','equilibrium_error','equilibrium',zstar);
sol.extdata=struct();S.solutions={sol};fprintf('LOCAL_ERROR_COORDINATE_DONE wall=%.2f steps=%d\n',S.wallSeconds,numel(sol.x)-1);
 function dq=rhs(~,q)
  X=reshape(q(P.iX),10,5).';V=reshape(q(P.iV),2,5).';[A,b,un,aux]=audit_core('rows',Pe,X,V,true);
  w=q([P.iU P.iZ]);la=q(P.iL);g=A*w+b;dl=g;dl(la<=0)=max(g(la<=0),0);
  dw=-H*w+[un;zeros(25,1)]-A.'*la;acc=aux.f0+reshape(q(P.iU),2,5).';
  dq=[reshape(aux.Xdot.',50,1);reshape(acc.',10,1);dw/P.epsilon;dl/P.epsilon];
 end
 function J=jac(t,q)
  J=zeros(P.dim);X=reshape(q(P.iX),10,5).';V=reshape(q(P.iV),2,5).';[A,b]=audit_core('rows',Pe,X,V,true);active=q(P.iL)>0 | A*q([P.iU P.iZ])+b>0;
  for j=1:60,h=1e-7;qp=q;qm=q;qp(j)=qp(j)+h;qm(j)=qm(j)-h;J(:,j)=(rhs(t,qp)-rhs(t,qm))/(2*h);end
  J(P.iV,P.iU)=eye(10);J([P.iU P.iZ],[P.iU P.iZ])=-H/P.epsilon;J([P.iU P.iZ],P.iL)=-A.'/P.epsilon;J(P.iL,[P.iU P.iZ])=(A.*active)/P.epsilon;
 end
end
