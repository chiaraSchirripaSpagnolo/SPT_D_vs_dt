function [nn,ctrs,edges,h] = hist3werr(varargin)
%HIST3WERR Three-dimensional histogram of bivariate data with weights and uncertainties.
%   [nn,ctrs,edges,h] = HIST3WERR(varargin)
%
%   HIST3WERR(X) bins the elements of the columns 1 and 3 of the  M-by-5 matrix
%   X into a 10-by-10 grid of equally-spaced containers. 
%   Each column of X corresponds to one dimension in the bin grid.
%   
%   X = [x1, sigma1, x2, sigma2, weights]
%
%   ... = HIST3WERR(...,'makeplot',1,...): a 2D-histogram is plotted.
%
%   HIST3WERR(X,NBINS) plots a histogram using an NBINS(1)-by-NBINS(2) grid of
%   bins.  HIST3WERR(X,'Nbins',NBINS) is equivalent to HIST3WERR(X,NBINS).
%
%   HIST3WERR(X,CTRS), where CTRS is a two-element cell array of numeric
%   vectors with monotonically non-decreasing values, uses a 2D grid of
%   bins centered on CTRS{1} in the first dimension and on CTRS{2} in the
%   second.  HIST3WERR assigns rows of X falling outside the range of that grid
%   to the bins along the outer edges of the grid, and ignores rows of X
%   containing NaNs.  HIST3WERR(X,'Ctrs',CTRS) is equivalent to HIST3WERR(X,CTRS).
%
%   HIST3WERR(X,'Edges',EDGES), where EDGES is a two-element cell array
%   of numeric vectors with monotonically non-decreasing values, uses a 2D
%   grid of bins with edges at EDGES{1} in the first dimension and at
%   EDGES{2} in the second.  The (i,j)-th bin includes the value X(k,:) if
%
%      EDGES{1}(i) <= X(k,1) < EDGES{1}(i+1) and
%      EDGES{2}(j) <= X(k,2) < EDGES{2}(j+1).
%
%   Rows of X that that fall on the upper edges of the grid, EDGES{1}(end)
%   or EDGES{2}(end), are counted in the (I,j)-th or (i,J)-th bins, where
%   I and J are the lengths of EDGES{1} and EDGES{2}.  HIST3WERR does not count
%   rows of X falling outside the range of the grid.  Use -Inf and Inf in
%   EDGES to include all non-NaN values.
%
%   N = HIST3WERR(X,...) returns a matrix containing the number of elements of
%   X that fall in each bin of the grid, and does not plot the histogram.
%   
%   [N,C] = HIST3WERR(X,...) returns the positions of the bin centers in a
%   1-by-2 cell array of numeric vectors, and does not plot the histogram.
%
%   HIST3WERR(AX,X,...) plots into AX instead of GCA.
%
%   HIST3WERR(..., 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   graphics parameter name/value pairs to fine-tune the plot.
%
%   see also HIST3

%   Copyright 2011 Luin; adapted from hist3 (Copyright 1993-2009 The
%   MathWorks, Inc.)
%   $Revision: 0.1 $  $Date: 2011/10/14 16:01:00 $


%nargs=length(varargin);
%args=varargin;
[cax,args,nargs] = axescheck(varargin{:});

if nargs < 1
    error('stats:hist3werr:TooFewInputs', 'Requires X.')
end
xdata = args{1};
ndata = size(xdata,1);
switch size(xdata,2)
    case 5
    case 0
    case 1
        xdata(:,[2,3,4,5])=repmat([0,0,0,1],ndata,1);
    case 2
        xdata=[xdata(:,1),zeros(ndata,1),xdata(:,2),zeros(ndata,1),ones(ndata,1)];
    case 3
        xdata=[xdata(:,1),zeros(ndata,1),xdata(:,2),zeros(ndata,1),xdata(:,3)];
    case 4
        xdata(:,5)=ones(ndata,1);
end;

    
xdata=xdata((~any(isnan(xdata),2)&all(isfinite(xdata(:,[2,4,5])),2)),:);
x=xdata(:,[1,3]);
weights=real(xdata(:,5));
sigmas=real(xdata(:,[2,4]));
sigmanonzero=any(sigmas,1);
% See if nbins/ctrs was given as the second argument, or only name/value
% pairs.
if nargs > 1 && ~ischar(args{2})
    binSpec = args{2};
    args = args(3:end); % strip off x and nbins/ctrs
else
    binSpec = [];
    args = args(2:end); % strip off x
end

% Process input parameter name/value pairs
pnames = {'nbins','ctrs','edges','makeplot','halfhalf','log1','log2'};
dflts =  { [],     [],       [],    0,         0,         0,     0};
if verLessThan('matlab', '8')
    [errid,errmsg,nbins,ctrs,edges,makeplot,half,log1,log2,plotArgs] = internal.stats.getargs(pnames, dflts, args{:});
    if ~isempty(errmsg)
        error(['Luin:histwerr2D:' errid], errmsg);
    end
else
    [nbins,ctrs,edges,makeplot,half,log1,log2,~,plotArgs] = internal.stats.parseArgs(pnames, dflts, args{:});
end
%{
if log1, fun{1}=@(x,mu,s) normcdfc(real(10.^x),real(10.^mu),sigma);
else fun{1}=@(x,mu,s) normcdfc(x,mu,sigma); end;

if log2, fun{2}=@(x,mu,s) normcdfc(real(10.^x),real(10.^mu),sigma);
else fun{2}=@(x,mu,s) normcdfc(x,mu,sigma); end;
%}
if log1, fun{1}=@(x) real(10.^x);
else fun{1}=@(x) real(x); end;

if log2, fun{2}=@(x) real(10.^x);
else fun{2}=@(x) real(x); end;

% Make sure they haven't mixed 'nbins'/'ctrs'/'edges' name/value pairs with
% the CTRS or NBINS unnamed second arg syntax, or used more than one of
% those parameter name/value pairs.
if (isempty(nbins)+isempty(ctrs)+isempty(edges)+isempty(binSpec)) < 3
    error('stats:hist3werr:AmbiguousBinSpec', 'Ambiguous specification of bins.');
elseif ~isempty(binSpec)
    if iscell(binSpec)  % hist3werr(x,ctrs)
        ctrs = binSpec;
    else                % hist3werr(x,nbins)
        nbins = binSpec;
    end
end

if ~isempty(nbins)
    % Use the specified number of bars in each direction, centers and edges
    % to be determined.
    histBehavior = true;
    if ~(isnumeric(nbins) && numel(nbins)==2)
        error('stats:hist3werr:BadNbins', ...
              'The number of bins must be specified with a 2-element numeric vector.');
    end
    autobins = true;
    
elseif ~isempty(ctrs)
    % Use the specified bin centers.
    histBehavior = true;
    if ~(iscell(ctrs) && numel(ctrs)==2 && isnumeric(ctrs{1}) && isnumeric(ctrs{2}))
        error('stats:hist3werr:BadCtrs', ...
              'Bin centers must be specified with a cell array containing two numeric vectors.');
    end
    ctrs = {ctrs{1}(:)' ctrs{2}(:)'};
    autobins = false;
    nbins = [length(ctrs{1}) length(ctrs{2})];
    
elseif ~isempty(edges)
    % Use the specified bin edges.
    histBehavior = false;
    if ~(iscell(edges) && numel(edges)==2 && isnumeric(edges{1}) && isnumeric(edges{2}))
        error('stats:hist3werr:BadEdges', ...
              'Bin edges must be specified with a cell array containing two numeric vectors.');
    end
    edges = {edges{1}(:)' edges{2}(:)'};
    autobins = false;
    % Just as with histc, there will be #edges bins
    nbins = [length(edges{1}) length(edges{2})];
    
else
    % Assume a 10x10 grid of bars, centers and edges to be determined.
    histBehavior = true;
    autobins = true;
    nbins = [10 10];
end

[nrows,ncols] = size(x);
if ncols ~= 2
    error('stats:hist3werr:WrongNumCols', 'X must be a matrix with two columns.');
end

% Special case for empty data (follows what HIST does).
if isempty(x)
    if autobins
       ctrs = {1:nbins(1) 1:nbins(2)};
    end
    nn = zeros(nbins); % Nothing to count, return nbins(1) by nbins(2) zeros
    
else
    % Bin each observation in the x-direction, and in the y-direction.
    if ~all(sigmanonzero) && ~half
        bin = zeros(nrows,2);
    end;
    histcEdges=cell(1,2);
    binwidth=cell(1,2);
    for i = 1:2
        minx = min(x(:,i));
        maxx = max(x(:,i));
        
        % If only the number of bins was given, compute edges and centers
        % for equal-sized bins spanning the data.
        if autobins
            if isinf(minx) || isinf(maxx)
                error('stats:hist3werr:InfData', ...
                      'Bin centers or edges must be specified when data contain infinite values.');
            elseif minx == maxx
                minx = minx - floor(nbins(i)/2) - 0.5;
                maxx = maxx + ceil(nbins(i)/2) - 0.5;
            end
            binwidth{i} = (maxx - minx) / nbins(i);
            edges{i} = minx + binwidth{i}*(0:nbins(i));
            ctrs{i} = edges{i}(1:nbins(i)) + binwidth{i}/2;
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            histcEdges{i} = [-Inf edges{i}(2:end-1) Inf];
            
        % If the bin centers were given, compute their edges and widths.
        elseif histBehavior
            c = ctrs{i};
            dc = diff(c);
            edges{i} = [c(1) c] + [-dc(1) dc dc(end)]/2;
            binwidth{i} = diff(edges{i});
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            histcEdges{i} = [-Inf edges{i}(2:end-1) Inf];
            
        % If the bin edges were given, compute their widths and centers (if
        % asked for).
        else % if ~histBehavior
            e = edges{i};
            de = diff(e);
            histcEdges{i} = e;
            % Make the display mimic bar's histc behavior: an implied bin
            % above edges(end), the same width as the last explicit one.
            % ctrs, edges, and binwidth need that explicitly, histcEdges
            % doesn't.
            edges{i} = [e e(end)+de(end)];
            binwidth{i} = [de de(end)];
            if nargout > 1
                c = zeros(size(de));
                c(1) = e(1) + de(1)/2;
                for j = 2:length(c)
                    c(j) = 2*e(j) - c(j-1);
                end
                % When edges are specified, it may not be possible to return
                % centers for which the edges are midpoints.  Warn if that's
                % the case.
                if any(c <= e(1:end-1)) || ...
                   abs(c(end) - (e(end)-de(end)/2)) > 1000*eps(de(end));
                    warning('stats:hist3werr:InconsistentEdges', ...
                            'Cannot compute centers that are consistent with EDGES.');
                    c = e(1:end-1) + de/2;
                end
                ctrs{i} = [c e(end)+de(end)/2];
            end
        end
        
        % Get the 1D bin numbers for this column of x.  Make sure +Inf
        % goes into the nth bin, not the (n+1)th.
        if ~sigmanonzero(i) && ~half
            [~,bin(:,i)] = histc(fun{i}(x(:,i)),[-Inf,fun{i}(histcEdges{i}(2:end))],1);
            bin(:,i) = min(bin(:,i),nbins(i));
        end
    end
    
    % Combine the two vectors of 1D bin counts into a grid of 2D bin
    % counts.
    if all(~sigmanonzero) && ~half
        nn = accumarray(bin(all(bin>0,2),:),weights(all(bin>0,2)),nbins);
    else
%        tic
        ncum=squeeze(sum(bsxfun(@times, ...
            repmat(weights,1,nbins(1)).*normcdfc(repmat(fun{1}(histcEdges{1}(2:end)),nrows,1),repmat(fun{1}(x(:,1)),1,nbins(1)),repmat(sigmas(:,1),1,nbins(1))),...
            permute(normcdfc(repmat(fun{2}(histcEdges{2}(2:end)),nrows,1),repmat(fun{2}(x(:,2)),1,nbins(2)),repmat(sigmas(:,2),1,nbins(2))),[1,3,2])),1));
%        toc
%
        nsum=[ncum(1,:); diff(ncum,1,1)];
        nn=[nsum(:,1),diff(nsum,1,2)];
    end
end
%%
if makeplot
del = .001; % space between bars, relative to bar size
n = nn;

% Build x-coords for the eight corners of each bar.
xx = edges{1};
xx = [xx(1:nbins(1))+del*binwidth{1}; xx(2:nbins(1)+1)-del*binwidth{1}];
xx = [reshape(repmat(xx(:)',2,1),4,nbins(1)); NaN(1,nbins(1))];
xx = [repmat(xx(:),1,4) NaN(5*nbins(1),1)];
xx = repmat(xx,1,nbins(2));

% Build y-coords for the eight corners of each bar.
yy = edges{2};
yy = [yy(1:nbins(2))+del*binwidth{2}; yy(2:nbins(2)+1)-del*binwidth{2}];
yy = [reshape(repmat(yy(:)',2,1),4,nbins(2)); NaN(1,nbins(2))];
yy = [repmat(yy(:),1,4) NaN(5*nbins(2),1)];
yy = repmat(yy',nbins(1),1);

% Build z-coords for the eight corners of each bar.
zz = zeros(5*nbins(1), 5*nbins(2));
zz(5*(1:nbins(1))-3, 5*(1:nbins(2))-3) = n;
zz(5*(1:nbins(1))-3, 5*(1:nbins(2))-2) = n;
zz(5*(1:nbins(1))-2, 5*(1:nbins(2))-3) = n;
zz(5*(1:nbins(1))-2, 5*(1:nbins(2))-2) = n;

cax = newplot(cax);
holdState = ishold(cax);

% Plot the bars in a light steel blue.
cc = repmat(cat(3,.75,.85,.95), [size(zz) 1]);

% Plot the surface, using any specified graphics properties to override
% defaults.
h = surf(cax, xx, yy, zz, cc, 'tag','hist3', plotArgs{:});

if ~holdState
    % Set ticks for each bar if fewer than 16 and the centers/edges are
    % integers.  Otherwise, leave the default ticks alone.
    if (nbins(1)<16)
        if histBehavior && all(floor(ctrs{1})==ctrs{1})
            set(cax,'xtick',ctrs{1});
        elseif ~histBehavior && all(floor(edges{1})==edges{1})
            set(cax,'xtick',edges{1});
        end
    end
    if (nbins(2)<16)
        if histBehavior && all(floor(ctrs{2})==ctrs{2})
            set(cax,'ytick',ctrs{2});
        elseif ~histBehavior && all(floor(edges{2})==edges{2})
            set(cax,'ytick',edges{2});
        end
    end
    
    % Set the axis limits to have some space at the edges of the bars.
    dx = range(edges{1})*.05;
    dy = range(edges{2})*.05;
    set(cax,'xlim',[edges{1}(1)-dx edges{1}(end)+dx]);
    set(cax,'ylim',[edges{2}(1)-dy edges{2}(end)+dy]);
    
    view(cax,3);
    grid(cax,'on');
    set(get(cax,'parent'),'renderer','zbuffer');
end
end

for i=1:2
    ctrs{i}=ctrs{i}.';
    edges{i}=edges{i}.';
end

end

function res=normcdfc(x,mu,sigma)
res=normcdf(x,mu,sigma);
if any(isnan(res(:)))
    if any(sigma==0)
        res(sigma==0 & x==mu)=0.5;
        res(sigma==0 & x>mu)=1;
        res(sigma==0 & x<mu)=0;
    end;
    res(x==+Inf)=1;
    res(x==-Inf)=0;
    res(sigma==Inf)=0;
end;
return;
end
    
    
        
