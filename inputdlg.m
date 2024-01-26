function Answer=inputdlg(varargin)
persistent shipped % stores a handle to the shipped MATLAB funtion of the same name
persistent MaxQuest
% tic;
if isempty(shipped) || isempty(MaxQuest)||MaxQuest<2
    this = ['inputdlg' '.m']; % the name of function in MATLAB we are shadowing
    list = which(this, '-all'); % find all the functions which shadow it
    f = strncmp(list, matlabroot, length(matlabroot)); % locate 1st in list under matlabroot
    list = list{find(f, 1)}; % extract from list the exact function we want to be able to call
    here = cd(list(1:end-length(this))); % temporarily switch to the containing folder
    shipped = str2func(this(1:end-2)); % grab a handle to the function
    cd(here); % go back to where we came from
    prevun=get(0,'Units');
    set(0,'Units','pixels');
    scrhgt=get(0,'ScreenSize');
    set(0,'Units',prevun);
    clear prevun;
    MaxQuest=floor((scrhgt(4)-120)/47);
%     MaxQuest=15;
end
% toc;
drawnow; pause(0.02);
if nargin<1
  varargin{1}=getString(message('MATLAB:uistring:popupdialogs:InputDlgInput'));
end
if ~iscell(varargin{1})
    varargin{1}=cellstr(varargin{1});
end
NumQuest=length(varargin{1});
if NumQuest<=MaxQuest
    Answer=shipped(varargin{:});
    drawnow; pause(0.05);
else
    Answer=cell(1,NumQuest);
    %    FigHeight=30+47*NumQuest; %(NumQuest+2)*5+  23 +NumQuest*19 +  NumQuest*23
    if length(varargin)<2
        root={''};
    else
        root=varargin{2};
    end
    mic=ceil(NumQuest/MaxQuest);
    for ic=1:mic
        varargin2={};
        range=(ic-1)*MaxQuest+1:min(ic*MaxQuest,NumQuest);
        varargin2{1}=varargin{1}(range);
        varargin2(2)={[int2str(ic) '/' int2str(mic) ' ' root]};
        if length(varargin)>2
            if length(varargin{3})==1
                varargin2=[varargin2,varargin(3)]; %#ok<*AGROW>
            else
                varargin2=[varargin2,{varargin{3}(range)}];
            end
        end
        if length(varargin)>3
                varargin2=[varargin2,{varargin{4}(range)}];
        end
        if length(varargin)>4
                varargin2=[varargin2,varargin(5)];
        end
        Answer(range) = shipped(varargin2{:});
        drawnow; pause(0.05);
    end
end
