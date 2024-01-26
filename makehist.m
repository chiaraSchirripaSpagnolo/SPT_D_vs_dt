function [h,hlog]=makehist(t,choice,fname_,name_,answ)
%[h,hlog]=makehist(t,choice,fname_,name_)
%choice=input('D12: 0; L with fixed bins: 1; free bins: 2; choice? (0,1,2) ');
%    fname_=input('File name root? ','s');
%    name_=input('Physical Quantity? ','s');
% answ~=0, save to file; answ=2, make plot

if nargin<5, answ=1; end

if ~isscalar(choice) || iscell(choice)
    [h,hlog]=histwerr(t,choice);
else
    switch choice
        case 0
            [h,hlog]=histwerr(t,0,0,0,-6,0.2);
        case 1
            [h,hlog]=histwerr(t,0,20,0.5);
        case 2
            [h,hlog]=histwerr(t);
        otherwise
            [h,hlog]=histwerr(t,choice);
    end
            
end

if answ
    savehist(name_,fname_,hlog,h);
end;

if answ==2
    figure;
    semilogx(hlog(:,2),hlog(:,5));
end