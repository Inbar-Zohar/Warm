function y = decimate_array(x,r,dim,varargin)

if ~exist('dim','var') || isempty(dim); dim = 2; end
if dim == 1; x = x'; end

l = ceil(size(x,2)/r);
y = nan(size(x,1),l);
for ind=1:size(x,1)
    y(ind,:) = decimate(x(ind,:),r);
end
if dim == 1; y = y'; end
