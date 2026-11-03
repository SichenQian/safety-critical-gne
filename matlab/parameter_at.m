function v=parameter_at(p,index)
% The manuscript allows row-dependent tightening parameters.
% Scalar configurations reproduce the original implementation exactly.
if isscalar(p),v=p;else,v=p(index);end
end
