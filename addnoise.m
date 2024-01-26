function im_out=addnoise(im_out,gaussnoise,nwh,nbl,poiss,whl,bll,maxn,minn)
%im_out=ADDNOISE(im_out,gaussnoise,nwh,nbl,poiss,whl,bll,maxn,minn)
%gaussnoise: add a noise with gaussian distribution with sigma=gaussnoise
%nwh: "salt" noise, add whl in nwh pixels, 
%nbl: "pepper" noise, multiply by bll nbl pixels
%poiss: if different than 0, add poissonian noise per pixel: each pixel
%  value is substituted with a number calculated by multiplying by poiss a
%  random number with poissonian distribution having as average;
%  (pixel value)/poiss; if poiss is >1, this is then substituted with a
%  number having uniform distribution within an interval of width poiss
%  around the result.
%if maxn is different than 0, and is a natural number, all pixel values are
%  rounded to the nearest integer.
%If maxn is different than zero, pixel values are clipped at maxn
%minn: minimum value acceptable for pixel value.
% default:
%     gaussnoise=0;
%     nwh=0;
%     nbl=0;
%     poiss=0;
%     whl=500;
%     bll=0;
%     maxn=0;
%     minn=-Inf;

if nargin<2 || isempty(gaussnoise)
    gaussnoise=0;
end
if nargin<3 || isempty(nwh)
    nwh=0;
end
if nargin<4 || isempty(nbl)
    nbl=0;
end
if nargin<5 || isempty(poiss)
    poiss=0;
end
if nargin<6 || isempty(whl)
    whl=500;
end
if nargin<7 || isempty(bll)
    bll=0;
end
if nargin<8 || isempty(maxn)
    maxn=0;
end
if nargin<8 || isempty(minn)
    minn=-Inf;
end

szs=size(im_out);
im_out=im_out+gaussnoise*randn(szs);
whn=randi(numel(im_out),1,nwh);
bln=randi(numel(im_out),1,nbl);
im_out(bln)=im_out(bln)*bll;
im_out(whn)=im_out(whn)+whl;

if poiss
    im_out(im_out>0)=poiss*poissrnd(im_out(im_out>0)/poiss);
    if poiss>1
        im_out=im_out+round(-poiss/2+poiss*rand(szs));
    end
end

if maxn
    im_out(im_out>maxn)=maxn;
    if round(maxn)==maxn
        im_out=round(im_out);
    end
end
im_out(im_out<minn)=minn;
    

