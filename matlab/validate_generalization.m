function validate_generalization()
P=audit_core('parameters');Q=frozen_core('parameters');rng(281);err=0;
for k=1:20
 X=repmat(reshape(P.y0.',1,10),5,1)+.02*randn(5,10);V=.02*randn(5,2);
 [A,b,u]=audit_core('rows',P,X,V,true);[B,c,v]=frozen_core('rows',Q,X,V,true);
 err=max(err,max(abs([A(:)-B(:);b-c;u-v])));
end
assert(err==0);P.obs=P.obs(1:2,:);P.rObs=P.rObs(1:2);P=audit_core('configure',P);
[A,b,u]=audit_core('rows',P,X,V,true);
keep=[1:22 24:25 27:28 30:31 33:34 36:45];
assert(isequal(A,B(keep,:))&&isequal(b,c(keep))&&isequal(u,v));
root=paper_data_root();write_json(fullfile(root,'tables','generalization_validation.json'),struct('baselineRowMaxError',err,'twoObstacleRowSelectionExact',true,'newStateDimension',P.dim));
fprintf('VARIABLE_OBSTACLE_GENERALIZATION_PASS\n');
end
