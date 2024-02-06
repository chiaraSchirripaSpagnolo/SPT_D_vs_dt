function [r,rlog,s2r,s2rlog]=histwerr(xdata,dmin,dmax,dx,logmin,dlogx)
% HISTWERR
% [r,rlog]=histwerr(xdata,dmin,dmax,dx,logmin,dlogx)
% gives distribution of xdata(:,1) with std xdata(:,2) and weight
% xdata(:,3) from dmin to dmax with spacing dx
% r: matrix with columns: bin center, bin end, hist, cumulative distr.
% rlog: matrix with column: bin center exp(10), bin center, bin end exp(10),...
%       bin end, hist, cumulative distr.
% s2r,s2rlog: matrices with columns: variance and uncertainties of histogram bins,
%             uncertainties of normalized histogram bins.
%
% More details on the uncertainties in “Molecular insight on the altered membrane
% trafficking of TrkA kinase dead mutants”. R. Amodeo, R. Nifosì, C. Giacomelli,
% C. Ravelli, L. La Rosa, A. Callegari, M.L. Trincavelli, S. Mitola, S. Luin,
% and L. Marchetti, Biochimica et Biophysica Acta (BBA) - Molecular Cell Research
% 1867(2), 118614 (2019). Doi: 10.1016/j.bbamcr.2019.118614.

no_data=false;
if nargin==0 || isempty(xdata)
    nodata;
    %     return
end

sigma0ist=1;
islog=0;
isbincenter=0;
if iscell(xdata)
    sigma0ist=0;
    xdata=xdata{1};
    if isempty(xdata),xdata=0; end
end
if nargin>1 && ~isempty(dmin) && iscell(dmin)
    dmin=dmin{1};
    islog=1;
end
ndata=size(xdata,1);
    function nodata
        r=zeros(1,5);
        s2r=Inf(1,5);
        rlog=zeros(1,7);
        s2rlog=Inf(1,7);
        no_data=true;
        if ~exist('xdata','var')
            xdata=[];
        end
    end
if ~ndata
    nodata;
    %     return;
end

if nargin>1 && ~isempty(dmin), nnbins=dmin; end
if nargin==2 && length(dmin)>1
    isbincenter=1;
    if any(diff(dmin)<=0)
        error('LUIN:histwerr','In histwerr(xdata,ctrs), ctrs must be a scalar, a vector with strictly increasing values, or a cell containing one of them');
    end
    szd=size(xdata,2);
    if szd<3
        xdata(:,3)=1;
        if szd<2
            xdata(:,2)=0;
        end
    end
    isok= ~isnan(xdata(:,1)) & ~isnan(xdata(:,2)) & ~isnan(xdata(:,3));
    xdata=xdata(isok,:);
    ndata=size(xdata,1);
    if ~ndata
        nodata;
        %         return;
    end
    bincenters=dmin(1:end);
    binlimits=mean([dmin(1:end-1);dmin(2:end)]);
    %binlimits=[binlimits,2*bincenters(end)-binlimits(end)];
    %dx=dmin(end)-dmin(end-1);
    dxmax=2*(bincenters(end)-binlimits(end));
    dxmin=2*(binlimits(1)-bincenters(1));
    dx=mean([dxmin,dxmax]);
    dmax=dmin(end-1);
    dmin=dmin(1);
    
else %end of histwerr(x,bincenters)
    
    if (nargin<3 || isempty(dmin) || isempty(dmax) || dmin == dmax || ~isfinite(dmin) || ~isfinite(dmax))
        if ~no_data
            dmin=min(xdata(isfinite(xdata(:,1)),1));
            dmax=max(xdata(isfinite(xdata(:,1)),1));
        else
            dmin=0;
            dmax=0;
            dx=1;
        end
    end
    
    if dmax<dmin
        temp=dmax;
        dmax=dmin;
        dmin=temp;
        clear temp;
    end
    
    if (nargin<4 || isempty(dx) || dx == 0 || ~isfinite(dx))
        if nargin==2
            if nnbins>1, dx=(dmax-dmin)./(nnbins-1);
            else dx=2.001*(dmax-dmin); end
        elseif ~no_data
            dx=(dmax-dmin)./sqrt(ndata);
        else
            dx=1;
        end
    end
    
    szd=size(xdata,2);
    if szd<3
        xdata(:,3)=1;
        if szd<2
            if ~no_data
                if dx==0, dx=mean(xdata(isfinite(xdata(:,1)),1)); end
                if dx==0, dx=1; end
            else
                dx=1;
            end
            %        xdata(:,2)=min(abs(xdata(:,1)),dx)/10;
            xdata(:,2)=0;
        end
    end
    
    isok= isfinite(xdata(:,1)) & isfinite(xdata(:,2)) & isfinite(xdata(:,3));
    %xdata(isnotok,3)=0;
    xdata=xdata(isok,:);
    ndata=size(xdata,1);
    if ~ndata
        nodata;
        %         return;
    end
    if dx==0, dx=mean(xdata(:,2)); end
    if dx==0, dx=1; end
    
    dx=abs(dx);
    
    binlimits=(dmin-dx/2.0:dx:dmax+dx/2.0);
    bincenters=mean([[dmin-1.5*dx, binlimits];[binlimits,binlimits(end)+dx]]);
    
end %end of else histwerr(x,bincenters)

xdata(:,2)=abs(xdata(:,2));

if any(xdata(:,2)) || ~sigma0ist
    xdata((xdata(:,2)==0),2)=min(abs(xdata(xdata(:,2)==0,1))/10,dx/10);
    xdata(xdata(:,2)==0,2)=dx/100;
end

nbins=length(binlimits);

if ~no_data
    if all(xdata(:,2))
        binlimmat=repmat(binlimits,ndata,1);
        datamat=repmat(xdata(:,1),1,nbins);
        sigmamat=repmat(xdata(:,2),1,nbins);
        %dxred=(binlimmat-datamat)./sigmamat;
        if islog
            gaussdata=normcdf(real(10.^binlimmat),real(10.^datamat),sigmamat);
        else
            gaussdata = normcdf(binlimmat,datamat,sigmamat);
        end
        %gaussdata(dxred>4) = 1;
        %gaussdata(dxred<-4) = 0;
        clear binlimmat datamat sigmamat % dxred;
        
        gaussdata=bsxfun(@times,xdata(:,3),gaussdata);
        histcum=sum(gaussdata,1);
        gaussdata(:,end+1)=xdata(:,3);
        gausssquare=(gaussdata(:,1:end-1)-gaussdata(:,2:end)).^2;
        s2r=[sum(gaussdata(:,1).^2,1),sum(gausssquare,1)];
        %histcum=sum(bsxfun(@times,xdata(:,3),gaussdata),1);
        histcum(end+1)=sum(xdata(:,3));
        histbin = histcum-[0,histcum(1:end-1)];
        s2r=(s2r-histbin.^2./ndata).*ndata./(ndata-1);
    else
        histcum = arrayfun(@(blm)sum(xdata(xdata(:,1)<blm,3)),binlimits);
        s2r=arrayfun(@(blm,blm1)sum(xdata(xdata(:,1)>=blm & xdata(:,1)<blm1,3).^2),[-Inf,binlimits],[binlimits,+Inf]);
%         ndatas2r=arrayfun(@(blm,blm1)sum(xdata(:,1)>=blm & xdata(:,1)<blm1),[-Inf,binlimits],[binlimits,+Inf]);
        histcum(end+1)=sum(xdata(:,3));
        histbin = histcum-[0,histcum(1:end-1)];
%         s2r=(s2r-hist.^2./ndatas2r).*ndatas2r./(ndatas2r-1);
        s2r=(s2r-histbin.^2./ndata).*ndata./(ndata-1); %doing like this, it approaches poissonian variance for ndatas2r<<ndata, even if all weights are equal (and, possibly, =1)
        s2r(~isfinite(s2r))=histbin(~isfinite(s2r)).^2;
        s2r(s2r==0)=+Inf;
    end
else
    histcum=zeros(size(binlimits));
    histcum(end+1)=0;
end

s2r=[s2r',sqrt(abs(s2r))'];
if no_data
    if ~exist('histbin','var')
        histbin=[histcum(1),diff(histcum)]; %should be all 0...
    end
    histnorm=histbin;
    s2r(:,end+1)=0;
else
    if isbincenter %da fare una cosa del genere per il caso log.
        histnorm = histbin./(histcum(end)*[dxmin,diff(binlimits),dxmax]);
        s2r(:,end+1)=(s2r(:,2)./(histcum(end)*[dxmin,diff(binlimits),dxmax]'));
    else
        histnorm = histbin./(histcum(end)*dx);
        s2r(:,end+1)=(s2r(:,2)./(histcum(end)*dx));
    end
end
r=[bincenters',[binlimits,dmax+1.5*dx]',histbin',histcum',histnorm'];

if nargout>1
    if dmax>0
        maxlog=log10(dmax);
        if dmin>0
            minlog=log10(dmin);
        else
            minlog=log10(min(xdata(xdata(:,1)>0,1)));
            if isempty(minlog) || ~isreal(minlog) || ~isfinite(minlog)
                minlog=maxlog;
            end
        end
    else
        minlog=0; maxlog=0;
    end
    nologmin=false;
    if nargin<5 || isempty(logmin)
        nologmin=true;
        logmin=minlog;
    end
    if logmin>maxlog %something strange, but let's avoid errors!
        temp=logmin;
        logmin=maxlog;
        maxlog=temp;
        clear temp;
    end
    if nargin<6 || isempty(dlogx)
        if no_data
            dlogx=1;
        else
            dlogx=(maxlog-logmin)./sqrt(ndata);
        end
    end
    if dlogx==0 && ~no_data, dlogx=dx/mean(abs(xdata(:,1))); end
    if isempty(dlogx) || dlogx==0 || ~isfinite(dlogx), dlogx=dx/100; end
    if nologmin
        logmin=logmin-dlogx;
    end
    
    logbinlimitsexp=(logmin+0.5*dlogx:dlogx:maxlog+0.5*dlogx);
    logbincentersexp=mean([[logmin-0.5*dlogx, logbinlimitsexp];[logbinlimitsexp,logbinlimitsexp(end)+1.0.*dlogx]]);
    
    logbinlimits=real(10.^logbinlimitsexp);
    logbincenters=real(10.^logbincentersexp);
    
    nlogbins=length(logbinlimits);
    
    if ~no_data
        if all(xdata(:,2))
            binlimmat=repmat(logbinlimits,ndata,1);
            datamat=repmat(xdata(:,1),1,nlogbins);
            sigmamat=repmat(xdata(:,2),1,nlogbins);
            %    dxred=(binlimmat-datamat)./sigmamat;
            gausslogdata = normcdf(binlimmat,datamat,sigmamat);
            %    gausslogdata(dxred>4) = 1;
            %    gausslogdata(dxred<-4) = 0;
            clear binlimmat datamat sigmamat %dxred;
            %             histlogc=sum(bsxfun(@times,xdata(:,3),gausslogdata),1);
            %         else
            %             histlogc = arrayfun(@(blm)sum(xdata(xdata(:,1)<blm,3)),logbinlimits);
            %         end
            %%
            gausslogdata=bsxfun(@times,xdata(:,3),gausslogdata);
            histlogc=sum(gausslogdata,1);
            gausslogdata(:,end+1)=xdata(:,3);
            gausssquare=gausslogdata(:,1:end-1).^2+gausslogdata(:,2:end).^2-2*gausslogdata(:,1:end-1).*gausslogdata(:,2:end);
            s2rlog=[sum(gausslogdata(:,1).^2,1),sum(gausssquare,1)];
            %histlogc=sum(bsxfun(@times,xdata(:,3),gausslogdata),1);
        else
            histlogc = arrayfun(@(blm)sum(xdata(xdata(:,1)<blm,3)),logbinlimits);
            s2rlog=arrayfun(@(blm,blm1)sum(xdata(xdata(:,1)>=blm & xdata(:,1)<blm1,3).^2),[-Inf,logbinlimits],[logbinlimits,+Inf]);
        end
        histlogc(end+1)=histcum(end);
        histlog = histlogc-[0,histlogc(1:end-1)];
        histlognorm = histlog./(histlogc(end)*dlogx);
        
        histbin = histlogc-[0,histlogc(1:end-1)];
        s2rlog=s2rlog-histbin.^2/ndata;
        s2rlog=[s2rlog',sqrt(abs(s2rlog))',sqrt(abs(s2rlog))'./(histlogc(end)*dlogx)];
    else
        histlogc=zeros(size(logbinlimits));
        histlogc(end+1)=0;
        histlog=histlogc;
        histlognorm=histlog;
        s2rlog=[histlog'.^2,histlog',histlognorm'];
    end
    
    rlog=[logbincentersexp',logbincenters',[logbinlimitsexp,logbinlimitsexp(end)+1.0.*dlogx]',...
        [logbinlimits,real(10^(logbinlimitsexp(end)+1.0.*dlogx))]',histlog',histlogc',histlognorm'];
end
return
end
