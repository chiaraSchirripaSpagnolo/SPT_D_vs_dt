function out=pref(act,group,name,value,varargin)
persistent strct
mlock;
if isempty(strct), strct=struct(); end
act=lower(act);
if nargin>=3 && ischar(name) && name(end)=='.'
    isstrct=true;
    name(end)='';
else
    isstrct=false;
end 
if strcmp(act,'add') && isfield(strct,group)&&isfield(strct.(group),name)&&...
        (~isstrct||(isstrct && (nargin>4 && isfield(strct.(group).(name),value))))
    error('Luos:pref','A preference with that GROUP and NAME already exists')
end

switch act
    case 'is'
        try
            switch nargin
                case 1
                    out=~isempty(fields(strct));
                case 2
                    out=~isempty(strct)&&all(isfield(strct,group));
                case 3
                    out=~isempty(strct)&&isfield(strct,group)&& all(isfield(strct.(group),name));
                case 4
                    if ~isstrct
                        out=~isempty(strct)&&isfield(strct,group)&& ...
                            isfield(strct.(group),name) && isequal(strct.(group).(name),value);
                    else
                        if ~isstruct(value)
                            out=~isempty(strct)&&isfield(strct,group)&& ...
                                isfield(strct.(group),name) && all(isfield(strct.(group).(name),value));
                        else
                            out=~isempty(strct)&&isfield(strct,group)&& ...
                                isfield(strct.(group),name) && all(isfield(strct.(group).(name),fieldnames(value)));
                            if out
                                for fnm=fieldnames(value)'
                                    out=out&&isequal(strct.(group).(name).(fnm{1}),value.(fnm{1}));
                                end
                            end
                        end
                    end
                case 5
                    out=~isempty(strct)&&isfield(strct,group)&& isfield(strct.(group),name) ...
                        && isfield(strct.(group).(name),value) && isequal(strct.(group).(name).(value),varargin{1});
                    
            end
            
        catch
            out=false;
        end
    case {'add','set'}
        switch nargin
            case 4
                if ~isstrct ||  ~strcmp(act,'add') || ~pref('is',group,name) || ~isstruct(strct.(group).(name))
                    strct.(group).(name)=value;
                else
                    strct.(group).(name)=copystruct(strct.(group).(name),value);
                end
            case 5
                if isstrct
                    strct.(group).(name).(value)=varargin{1};
                end
        end
     case 'rm'
        switch nargin
            case 1
                clear strct;
                strct=struct();
            case 2
                strct=rmfield(strct,group);
            case 3
                strct.(group)=rmfield(strct.(group),name);
            case 4
                if isstrct
                    strct.(group).(name)=rmfield(strct.(group).(name),value);
                end
        end
    case 'get'
        switch nargin
            case 1
                out=strct;
            case 2
                out=strct.(group);
            case 3
                out=strct.(group).(name);
            case 4
                if isstrct
                    out=strct.(group).(name).(value);
                else
                    if isfield(strct,group)&&isfield(strct.(group),name)
                        out=strct.(group).(name);
                    else
                        out=value;
                        strct.(group).(name)=value;
                    end
                end
            case 5
                if isfield(strct,group)&&isfield(strct.(group),name)&&isfield(strct.(group).(name),value) 
                    out=strct.(group).(name).(value);
                else
                    out=varargin{1};
                    strct.(group).(name).(value)=out; %error if strct.(group).(name) is not a struct!
                end
        end
    case 'uiget'
        if isfield(strct,group)&&isfield(strct.(group),name)&& ...
                ~strcmp(strct.(group).(name),'ask')
            out=strct.(group).(name);
        else
            if ispref(group,name)
                temp=getpref(group,name);
                rmpref(group,name);
            end
            out=uigetpref(group,name,value,varargin{:});
            strct.(group).(name)=getpref(group,name);
            if exist('temp','var')
                setpref(group,name,temp);
            else
                rmpref(group,name);
            end
        end
end