function varargout=importdataL(varargin)
persistent shipped % stores a handle to the shipped MATLAB funtion of the same name

if isempty(shipped)
    this = ['importdata' '.m']; % the name of function in MATLAB we are shadowing
    list = which(this, '-all'); % find all the functions which shadow it
    f = strncmp(list, matlabroot, length(matlabroot)); % locate 1st in list under matlabroot
    list = list{find(f, 1)}; % extract from list the exact function we want to be able to call
    here = cd(list(1:end-length(this))); % temporarily switch to the containing folder
    shipped = str2func(this(1:end-2)); % grab a handle to the function
    cd(here); % go back to where we came from
end
try
    [t.data,~,t.textdata,delm]=readfile(varargin{1});
    varargout{1}=t;
    varargout{2}=delm;
    varargout{3}=size(t.textdata,1)-size(t.data,1);
catch ME %#ok<NASGU>
    [varargout{1:nargout}]=shipped(varargin{:});
end
