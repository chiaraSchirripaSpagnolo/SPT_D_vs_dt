function savehist(name,fname,datalog,data,fit)
% SAVEHIST
% savehist(name,fname,[datalog,data],[fit])
% write comma separated histograms of name on files
% fname_name_[log]hist.csv
% fit in the last column of 'data'
if pref('is','Luos','savehist')&&~pref('get','Luos','savehist')
    return;
end
if nargin<3
    disp('savehist called with too few arguments: savehist(name,fname,datalog[,data[,fit]])');
    return;
end

% r: bin center, bin end, hist, cumulative distr.
% rlog: bin center exp(10), bin center, bin end exp(10), bin end, hist, cumulative distr.

if ~isempty(datalog)
    header={[name '_logcenter'],[name '_center'],[name '_logend'],[name '_end'],...
        [name '_freq'],[name '_cumulative'],[name '_freq_norm']};
    if size(datalog,2)==10
        header=[header,{'s2_freq','D_freq','D_freq_norm'}];
        datalog=datalog(:,[1:5,9,8,6,7,10]);
        header=header([1:5,9,8,6,7,10]);
    end
    csvwriteh([fname '_' name '_loghist.csv'],header,datalog);
end

if nargin>=4&&~isempty(data)
    header={[name '_center'],[name '_end'],...
        [name '_freq'],[name '_cumulative'],[name '_freq_norm']};
    if size(data,2)==8
        header=[header,{'s2_freq','D_freq','D_freq_norm'}];
        data=data(:,[1:3,7,6,4,5,8]);
        header=header([1:3,7,6,4,5,8]);
    end
    if nargin>=5 &&~isempty(fit)&&size(fit,1)==size(data,1)
        header=[header,{'fit'}];
        data=[data,fit];
    end
    csvwriteh([fname '_' name '_hist.csv'],header,data);
end

return;


    