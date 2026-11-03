function result=RUN_FINAL_REVISION(mode)
% Compatibility entry point for the numerical-revision package.
if nargin<1,mode='figures';end
repo=fileparts(fileparts(mfilename('fullpath')));addpath(repo);
result=reproduce_paper(mode);
end
