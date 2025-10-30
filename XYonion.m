function split = XYonion(X, y, fraction, nonion, mahalanobis)
% Apply splitcaltest using the 'onion' algorithm
% INPUT:
%    X = array of data (samples in rows)
%    y = vector of corresponding reference values
%    fraction = fraction of samples in calibration
%    
%    optional inputs
%       nonion = number of onion layers (default = 3)
%       mahalanobis = calculate Mahalanobis distances (= 1) instead of 
%       Euclidean (default = 0)
%
% OUTPUT:
%    split: vector indicating which samples correspond to calibration set
%
% Author: Jokin Ezenarro
% v1.0 (10/6/2025)
% Ref.: https://doi.org/10.1016/j.aca.2025.344229

% OPTIONS
if nargin<5
    mahalanobis = 0; % Do not calculate Mahalanobis distances
end
if nargin<4
    nonion = 3; % Number of onion layers
end
if nargin<3
    fraction = 0.66; % Fraction of samples in calibration
end
if nargin <2
    error('Not enough input arguments')
end

[m,~]  = size(X);
i0 = [1:m]; % indices for samples
split  = zeros(1,m);  % start assigning all samples to Cal

if mahalanobis==1
% Mahalanobis
    % Compute the covariance matrix of the dataset
    cov_matrix = cov(X);

    % Compute eigen decomposition
    [V, D] = eig(cov_matrix);

    % Construct the whitening transformation
    cov_matrix_inv_sqrt = V * diag(1 ./ sqrt(diag(D))) * V';

    % Apply the transformation
    X = X * cov_matrix_inv_sqrt;
    % Now Euclidean distances in X_transformed correspond to Mahalanobis distances in X
end

%ncalloop  = round(m/nonion*fraction); 
%ntestloop = round(m/nonion*(1-fraction));
loopfraction = 0.1;

for j=1:nonion
  m0     = length(i0);  % number of samples not yet assigned to cal or test
  ncalloop  =  round(loopfraction*m0*fraction);
  ntestloop =  round(loopfraction*m0*(1-fraction));

  if m0>0
    isel = distslct_xyonion(X(i0,:),y(i0,:),min(m0, ncalloop));
    % isel = distslct(X(i0,:),min(m0, ncalloop));
    split(i0(isel)) = 1; %Cal
    i0(isel)=[];
  end
  m0     = length(i0);
  if m0>0
    isel = distslct_xyonion(X(i0,:),y(i0,:),min(m0, ntestloop));
    % isel = distslct(X(i0,:),min(m0, ncalloop));
    split(i0(isel)) = 0; %Test
    i0(isel)=[];
  end
end

% randomly assign remaining samples to cal and test in desired proportion
m0 = length(i0);
nc = ceil(m0*fraction);
nt = m0-nc;
irem=[repmat(1,nc,1); repmat(0,nt,1)];
irem=irem(randperm(length(irem)));
split(i0) = irem';
split = logical(split)';

    function isel = distslct_xyonion(x,y,nosamps)
    %DISTSLCT Selects outter samples of data space.
    %  Inputs are the M by N data matrix (X) (class "double"), the 
    %  corresponding vector of reference values (y) and the number of 
    %  samples to select (nosamps).
    %
    %  Output (isel) is a vector of row indices for the selected samples.

    [q,n]         = size(x);
    x = x-mean(x,1); %mean center
    y = y-mean(y,1); %mean center
    isel          = zeros(1,nosamps);
    
    %stdslct gets n samples
	    if nosamps<n
           isel = stdsslct_xyonion(x,y,nosamps);
        else % when more samples need to be selected than there are variables
            isel(1:n) = stdsslct_xyonion(x,y,n); %select as many samples as variables normally
            distX        = zeros(q,1);
            distY        = zeros(q,1);
            
            for ii=1:n
                xsel = x(isel(ii),:);
                ysel = y(isel(ii),1);
                
                distX      = distX + sqrt(sum((x-xsel(ones(q,1),:)).^2,2));
                distY      = distY + sqrt(sum((y-ysel(ones(q,1),1)).^2,2));
                
                distX(isel(1:ii),:) = zeros(ii,1);
                distY(isel(1:ii),:) = zeros(ii,1);
            end
            
            dist = distX./max(distX) + distY./max(distY); %calculate SPXY distance
            [~,isel(n+1)] = max(dist); %first sample is furthest from the mean
            
            xsel = x(isel(n+1),:);
            ysel = y(isel(n+1),1);
            
            distX(isel(n+1)) = 0;
            distY(isel(n+1)) = 0;
            
            for ii=n+2:nosamps
                distX = distX + sqrt(sum((x-xsel(ones(q,1),:)).^2,2));
                distY = distY + sqrt(sum((y-ysel(ones(q,1),1)).^2,2));
                
                dist = distX./max(distX) + distY./max(distY); %calculate SPXY distance
                [~,isel(ii)] = max(dist);
                
                xsel = x(isel(ii),:);
                ysel = y(isel(ii),:);
                
                distX(isel(1:ii)) = zeros(ii,1);
                distY(isel(1:ii)) = zeros(ii,1);
            end
	    end
    end

    function isel = stdsslct_xyonion(x,y,nosamps)

        subset  = zeros(1,nosamps);
        idx = (1:length(y))'; % map to original row indices
        for i=1:nosamps   
            distX = sum(x.^2,2); %calculate distance to zero for X (norm)
            distY = sum(y.^2,2); %calculate distance to zero for y (norm)
            dist = distX./max(distX) + distY./max(distY); %calculate SPXY distance
            
            [~, k] = max(dist); % k is index in reduced set
            subset(1,i) = idx(k); % store ORIGINAL index

            rx0 = x(k,:); %save the furthest sample
            ry0 = y(k,:);

            % remove selected row from working copies and index map
            x(k,:) = [];
            y(k,:) = [];
            idx(k) = [];

            %make samples orthogonal to selected sample
            x = x - ((rx0*x')./(rx0*rx0'))' * rx0;
            y = y - ry0;
        end    
        isel = subset;
    end
end