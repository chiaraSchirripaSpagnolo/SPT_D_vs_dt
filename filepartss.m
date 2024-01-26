function parts=filepartss(str)
if isempty(str)||~ischar(str)
    parts={''};
    return;
end
parts={'','',''};
[parts{:}]=fileparts(str);
if ~isempty(parts{1}) && (parts{1}(end)==filesep||parts{1}(end)=='\'||parts{1}(end)=='/')
    a=regexp(str,'(\\\\|//)');
    if ~isempty(a)
        a=a(end);
        if strcmp(parts{1},str(1:a))
            parts{1}(end+1)=parts{1}(end);
        end
    end
end
if ~isempty(parts{1})&&(~(isempty(parts{2})&&isempty(parts{3})))
    parts=[dirparts(parts{1}),{[parts{2:3}]}];
elseif ~isempty(parts{1})
    parts=dirparts(parts{1});
else
    parts={[parts{2:end}]};
end
end

function parts=dirparts(str)
if isempty(str)
    parts={}; %should not be reached...
    return;
end
parts={'',''};
[parts{:},ext]=fileparts(str);
if ~isempty(ext)
    parts{2}=[parts{2},ext];
end
if ~isempty(parts{1})&&(parts{1}(end)==filesep||parts{1}(end)=='\'||parts{1}(end)=='/')
    a=regexp(str,'(\\\\|//)');
    if ~isempty(a)
        a=a(end);
        if strcmp(parts{1},str(1:a))
            parts{1}(end+1)=parts{1}(end);
        end
    end
end

if ~(isempty(parts{1})||isempty(parts{2}))
        parts=[dirparts(parts{1}),parts(2)];
elseif isempty(parts{1})
    parts=parts{2};
else
    parts=parts{1};
end
end
