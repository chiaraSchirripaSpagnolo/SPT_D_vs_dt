function fid=openw(nameout,cond)
%OPENW open file with name nameout and rule cond.
%loop until it can be opened, eventually using the new file name entered.
    while(1)
        fid=fopen(nameout,cond);
        if fid~=-1, break; end
        if any(lower(cond)=='w'|lower(cond)=='a')
            ddir=fileparts(nameout);
            ddirs=filepartss(nameout);
            if isempty(ddirs{1}) || ~(ddirs{1}(end)==filesep || ddirs{1}(end)=='\' || ddirs{1}(end)=='/')
                ddir=fullfile(pwd,ddir);
            end
            if ~isempty(ddir)&& ~exist(ddir,'dir')
                mkdir(ddir);
                continue;
            end
        end
        s = input(['cannot open file ' regexptranslate('escape',nameout) ' in dir ' regexptranslate('escape',pwd) ' (? to stop)'],'s');
        if ~isempty(s)
            if s=='?', break; end
            nameout = s;
        end
    end
end