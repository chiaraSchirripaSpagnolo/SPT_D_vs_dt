function str=catstruct(str1,str2,dim)
if nargin<3||isempty(dim)
  dim = find(size(str2)~=1, 1 );
  if isempty(dim)
      dim = find(size(str1)~=1, 1 );
    if isempty(dim), dim = 1; end
  end
end
try
    str=cat(dim,str1,str2);
catch %#ok<CTCH>
    if isempty(str1)
        str=str2;
    elseif isempty(str2)
        str=str1;
    else
        str=str1;
        szs=length(str1);
        for f = fieldnames(str2)'
            str(szs+1).(f{1}) = str2(1).(f{1});
        end
        if szs==1 && dim==1
            str=str';
        end
        
        if length(str2)>1
            for f = fieldnames(str1)'
                str2(1).(f{1}) = [];
            end
            str=cat(dim,str,str2(2:end));
        end
    end
end
