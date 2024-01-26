function im_out=imagespots(szs,centers,widths,amp)
if nargin<1 || isempty(szs)
    szs=[100,100];
end
if nargin<2 || isempty(centers)
    centers=[50,50];
end
if nargin<3 || isempty(widths)
    widths=[3,3];
end
if nargin<4 || isempty(amp)
    amp=1;
end
if numel(szs)==1
    szs=[szs,szs];
end
szs=szs(1:2);
if numel(centers)<2
    warning('LUIN:imagespots','centers should be a 2-columns array');
    centers=[centers,centers];
end
if size(centers,2)<2
    centers=centers';
end
centers=centers(:,1:2);
if size(widths,2)==1
    if size(widths,1)==2 && size(centers,1)==1
        widths=widths';
    else
        widths=[widths,widths];
    end
end
widths=widths(:,1:2);
amp=amp(:);
sv=[size(centers,1),size(widths,1),size(amp,1)];
sv1=sv==1;
tot1=sum(sv1);
varvec={'centers','widths','amp'};
switch(tot1)
    case 0
        if any(diff(sv))
            error('LUIN:imagespots',['wrong length of input arrays ',int2str(sv)]);
        end
    case 1
        if diff(sv(~sv1))
            error('LUIN:imagespots',['wrong length of input arrays ',int2str(sv)]);
        else
            eval([varvec{sv1} '=repmat(' varvec{sv1} ',sv(find(~sv1,1)),1);']);
        end
    case 2
            eval([varvec{find(sv1,1,'first')} '=repmat(' varvec{find(sv1,1,'first')} ',sv(~sv1),1);']);
            eval([varvec{find(sv1,1,'last')} '=repmat(' varvec{find(sv1,1,'last')} ',sv(~sv1),1);']);
end
im_out=hist3werr([centers(:,1),widths(:,1),centers(:,2),widths(:,2),amp],{0:szs(1)+1,0:szs(2)+1});
im_out(1,:)=[];
im_out(end,:)=[];
im_out(:,1)=[];
im_out(:,end)=[];