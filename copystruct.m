function [sout,mfields] = copystruct(sin,s,onlyex)
%COPYSTRUCT copies the struct s in sin, overriding the fields in sin.

if nargin<2 || ~isequal(size(sin),size(s))
    error('Luin:copystruct','copystruct must be called with two struct arrays of the same size');
end
if nargin<3
    onlyex=0;
end
sout = sin;
mfields={};
for f = fieldnames(s)'
    try
        if ~isfield(sin,f)
            mfields=[mfields,f]; %#ok<AGROW>
            if onlyex
                continue;
            end
        end
        [sout.(f{1})] = s.(f{1});
    catch ME
        if ~(strcmp(ME.identifier,'MATLAB:class:SetProhibited') || ...
                strcmp(ME.identifier,'lccb:set:readonly'))
            rethrow(ME);
        end
    end
end