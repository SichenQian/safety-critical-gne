function validate_error_coordinates()
P=design_parameters();Pe=P;Pe.errorCoordinates=true;rng(32);err=zeros(1,4);
Xs=repmat(reshape(P.ystar.',1,10),5,1);
for k=1:30
 X=1e-4*randn(5,10);V=1e-4*randn(5,2);[A,b,u,aux]=audit_core('rows',P,Xs+X,V,true);[Ae,be,ue,xe]=audit_core('rows',Pe,X,V,true);
 err=max(err,[max(abs(A-Ae),[],'all'),max(abs(b-be)),norm(u-ue,inf),max(abs([aux.Xdot(:)-xe.Xdot(:);aux.f0(:)-xe.f0(:)]))]);
end
assert(max(err)<1e-10);root=paper_data_root();
write_json(fullfile(root,'tables','error_coordinate_regression.json'),struct('maxDifferences',err,'tolerance',1e-10,'analyticIdentities',{{'L(y) times consensus equilibrium = 0','B ystar = p'}},'onlineQPSubstitution',false));
fprintf('ERROR_COORDINATE_REGRESSION_PASS %.3g\n',max(err));
end
