function dz=nominal_hb_rhs(P,z)
%NOMINAL_HB_RHS The 60-state nominal_HB_stacked equations in the manuscript.
% X(i,:) is agent i's estimate of all five planar decisions. Its own
% decision derivative is V(i,:); the other estimates follow consensus.
X=reshape(z(1:50),10,5).';V=reshape(z(51:60),2,5).';
Y=zeros(5,2);for i=1:5,Y(i,:)=X(i,2*i-1:2*i);end
[~,L]=audit_core('graph',P,Y);rho=P.c*L*X;Xdot=-rho;
acc=-P.gamma*V;
for i=1:5
 ii=2*i-1:2*i;Xdot(i,ii)=V(i,:);
 localY=reshape(X(i,:),2,5).';
 acc(i,:)=acc(i,:)-rho(i,ii)-(P.B(i,:)*localY-P.p(i,:));
end
dz=[reshape(Xdot.',50,1);reshape(acc.',10,1)];
end
