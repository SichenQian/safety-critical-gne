function validate_adopted_equations()
P=design_parameters();rng(20261101);err=0;
for k=1:25
 X=repmat(reshape(P.y0.',1,10),5,1)+.02*randn(5,10);V=.02*randn(5,2);
 [A,b,u,D]=audit_core('rows',P,X,V,true);[Ae,be,ue,E]=frozen_core('rows',P,X,V,true);
 err=max(err,max(abs([A(:)-Ae(:);b-be;u-ue;D.Xdot(:)-E.Xdot(:);D.f0(:)-E.f0(:)])));
end
assert(err==0);root=paper_data_root();
write_json(fullfile(root,'tables','adopted_equation_regression.json'),struct('randomStates',25,'maximumDifference',err,'comparison','Current adopted parameters in original frozen equation implementation versus generalized implementation'));
fprintf('ADOPTED_EQUATIONS_EXACT_PASS\n');
end
