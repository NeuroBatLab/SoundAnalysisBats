function H = hopkins_stats (X,minx,maxx,m,dist,n_neighbor,flip)
% hopkins_stats calculates H for Hopkins statistics estimation
%% set default values
if nargin<2
    minx = min(X);
end
if nargin<3
    maxx = max(X);
end
if nargin<4
    m=round(0.8 * size(X,1)); % sample from random distribution 10% of total number of element in X
end
%enforce a minimum of samples
if m<max(0.1*size(X,1),10*size(X,2)*2)
    m=max(0.1*size(X,1),10*size(X,2)*2);
end
if nargin<5
    dist = 'uniform';
end
if nargin<6
    n_neighbor = 1;
end
if nargin<7
    flip = 0;
end
%number of datapoints
N=size(X,1);
%range
Range=maxx-minx;
%% Get ready m samples from random distribution
ujd=nan(m,size(X,2));
bool=0;
if strcmp(dist,'uniform')
    %generate m random sampling points following a uniform distribution on each
    %dimension of X
    bool=1;
    for ii=1:size(X,2)
        ujd(:,ii)=minx(ii)+Range(ii)*rand(m,1);
    end
elseif strcmp(dist,'normal')
    %generate m random sampling points, normal
    %distribution based on mean and std of data points
    bool=1;
    for ii=1:size(X,2)
        ujd(:,ii)=normrnd(nanmean(X(:,ii)),nanstd(X(:,ii)),m,1);
    end
elseif strcmp(dist,'uniform_convex_hull')
    %generate m random sampling points in delaunay
    %triangulation of data points
    bool=0;
    r=0; % counter for the number of samples (needs to get m)
    DT=delaunayTriangulation(X);
    for rep=1:m*10
        sp_r = nan(1,size(X,2)); % random sample from random distribution
        for ii=1:size(X,2)
            sp_r(ii)=minx(ii)+Range(ii)*rand(1,1);
        end
        % finding if the random sample sp_r belong to the convex hull DT
        t=tsearchn(X,DT.ConnectivityList,sp_r);
        if ~isnan(t)
            r=r+1;
            ujd(r,:)=sp_r;
            if r==m
                bool=1;
                break
            end
        end
    end
    if bool==0
        error('Delaunay values error')
    end
else
    error('Distribution invalid')
end
%% Select m random data points from X without re-sampling:
Index=nan(m,1);% contains the indices of the m samples in X
wjd=nan(m,size(X,2)); % contains the values of the m samples taken from X
ii=0;
while ii <= m
    rn=round(N*rand);
    if ~any(Index == rn)
        if ~rn==0
            ii=ii+1;
            Index(ii)=rn;
            wjd(ii,:)=X(rn,:);
        end
    end
end
% figure;subplot(2,1,1);plot(ujd(:,1),ujd(:,2),'.')
% subplot(2,1,2);plot(wjd(:,1),wjd(:,2),'.')
%% Calculate distances between points
if bool
    U =nan(m,1); % average distance between each sample from random distribution and observed data
    W =nan(m,1); % average distance between each sample from observed distribution and observed data
    % minimal distance between random sampling points and data points:
    for ii = 1:m
        [~,D] = knnsearch(ujd(ii,:),X,'K',n_neighbor+1,'NSMethod','exhaustive');
        U(ii)=nanmean(D);
    end
    %U = (min(dist2 (X, ujd), [], 1));
    % minimal distance between random data points and data points:
    for ii = 1:m
        [~,D] = knnsearch(wjd(ii,:),X,'K',n_neighbor+1,'NSMethod','exhaustive');
        W(ii)=nanmean(D);
    end
    % leaveIndex = removerows((1:N)',Index);
    % dataR = X(leaveIndex,:);
    % % minimal distance between residual data points and random data points:
    % W = (min(dist2 (dataR, wjd), [], 1));
    if flip
        H=nansum(W)/(nansum(U)+nansum(W));
    else
        H=nansum(U)/(nansum(U)+nansum(W));
    end
else
    H=NaN;
end
end
