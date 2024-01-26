function A=scanfirsts(fname,n,headers,delimiter)

if nargin<2, error('too few args in scanfirsts'); end
if nargin<3, headers=0; end
if nargin<4, delimiter='\t'; end

fid = openw(fname,'rt');
for i=1:headers
    textscan(fid, '%s %*[^\n]',1);
end

formatstring=repmat('%f ',1,n);

C=textscan(fid, [formatstring '%*[^\n]'],'delimiter',delimiter,...
                                        'CollectOutput', 1);
A=C{1};
fclose(fid);