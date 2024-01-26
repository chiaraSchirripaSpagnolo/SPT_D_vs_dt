function [im_out,centers,im_nonoise]=randspots(szs,n,w,amp,gaussnoise,nwh,nbl,bg,varargin)
% [im_out,centers]=randspots(szs,n,w,amp,gaussnoise,nwh,nbl,bg,rel_da,rel_dw,poiss,whl,bll,maxn,minn)
% n: if double vector, n(j) number of particles of corresponding w(j) and a(j); if cell vector,
% n{j} contains the coordinates in a sv x 2 matrix
%gaussnoise: add a noise with gaussian distribution with sigma=gaussnoise
%nwh: "salt" noise, add whl in nwh pixels, 
%nbl: "pepper" noise, multiply by bll nbl pixels
%poiss: if different than 0, add poissonian noise per pixel: each pixel value is divided by a random number with poissonian distribution with average pixel value / poiss, multiplied by poiss.
%if maxn is different than 0, and is a natural number, all pixel values are
%rounded to the nearest integer.
%If maxn is different than zero, pixel values are clipped at maxn
%minn: minimum value acceptable for pixel value.
%
if nargin<1 || isempty(szs)
    szs=[100,100];
end
if nargin<2 || isempty(n)
    n=10;
elseif iscell(n)
    if iscell(n{1})
        centers=n{1};
    else
        centers=n;
    end
    n=cellfun(@(c) size(c,1),centers);
    n(n==2)=cellfun(@(c) size(c,2),centers(n==2));
else
    centers=cell(length(n),1);
end
if nargin<3 || isempty(w)
    w=[1,1];
end
if nargin<4 || isempty(amp)
    amp=1;
end
if nargin<5 || isempty(gaussnoise)
    gaussnoise=0;
end
if nargin<6 || isempty(nwh)
    nwh=0;
end
if nargin<7 || isempty(nbl)
    nbl=0;
end
if nargin<8 || isempty(bg)
    bg=0;
end
if numel(szs)==1
    szs=[szs,szs];
end
szs=reshape(szs(1:2),1,2);
n=n(:);
if numel(w)==numel(n)
    w=[w(:),w(:)];
elseif size(w,2)==1
    if size(w,1)==2 && size(n,1)==1
        w=w';
    else
        w=[w,w];
    end
end
w=w(:,1:2);
amp=amp(:);
sv=[size(n,1),size(w,1),size(amp,1)];
sv1=sv==1;
tot1=sum(sv1);
varvec={'n','w','amp'};
switch(tot1)
    case 0
        if any(diff(sv))
            error('LUIN:randspots',['wrong length of input arrays ',int2str(sv)]);
        end
    case 1
        if diff(sv(~sv1))
            error('LUIN:randspots',['wrong length of input arrays ',int2str(sv)]);
        else
            eval([varvec{sv1} '=repmat(' varvec{sv1} ',sv(find(~sv1,1)));']);
        end
    case 2
            eval([varvec{find(sv1,1,'first')} '=repmat(' varvec{find(sv1,1,'first')} ',sv(~sv1));']);
            eval([varvec{find(sv1,1,'last')} '=repmat(' varvec{find(sv1,1,'last')} ',sv(~sv1));']);
end
im_out=bg*ones(szs);

for cc=1:length(n)
    if ~isempty(varargin)&&~isempty(varargin{1})
        rel_da=abs(1+varargin{1}(1)*randn(n(cc),1));
    else
        rel_da=1;
    end
    if numel(varargin)>1 && ~isempty(varargin{2})
        rel_dw=sqrt(1+(bsxfun(@times,varargin{2}(min(cc,numel(varargin{2}))),randn(n(cc),2))).^2);
%         rel_dw=abs(1+bsxfun(@times,varargin{2}(min(cc,numel(varargin{2}))),randn(n(cc),2)));
    else
        rel_dw=1;
    end
    wv=bsxfun(@times,rel_dw,w(cc,1:2));
    rng('shuffle');
    if isempty(centers{cc})
        centers{cc}=bsxfun(@times,rand(n(cc),2),szs-1)+1;
    elseif size(centers{cc},2)~=2 && size(centers{cc},1)==2
        centers{cc}=centers{cc}';
    end
    im_out=im_out+imagespots(szs,centers{cc},wv,amp(cc)*2*pi*w(cc,1)*w(cc,2).*rel_da);
%     im_out=im_out+imagespots(szs,centers{cc},wv,amp(cc)*2*pi*wv(:,1).*wv(:,2).*rel_da);
end
if nargout>2
    im_nonoise=im_out;
end
if any([gaussnoise,nwh,nbl,[varargin{:}]])
    im_out=addnoise(im_out,gaussnoise,nwh,nbl,varargin{3:end});
end