function root=paper_data_root(value)
% Select the output tree without changing any numerical equation.
persistent active
if nargin,active=char(value);end
if isempty(active),root=fileparts(fileparts(mfilename('fullpath')));else,root=active;end
end
