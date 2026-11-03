function P=design_parameters()
% Contrast-enhanced candidate 15. Fixed geometry/game and complete online ODE.
P=audit_core('parameters'); % All three original obstacles and radii retained.
P.dc=2.2;P.dsafe=.345;P.lambdaLower=.04;P.c=200;P.gamma=65;
P.K1=4;P.K2=4.6;
P.DeltaInt=.5;P.DeltaObs=.7;P.DeltaC=.4;P.epsilon=5e-6;P.xi=1.4;
P.trackingBudget=.061;P.connMismatchBudget=.03;P=audit_core('configure',P);
end
