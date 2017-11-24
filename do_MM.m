function [ res_MM,minBIC ] = do_MM( data, fdrV, varargin )
%DO_MM Calculates the sparsified partial correlation matrix using a
%mixture modeling based approach from BOLD time series input datasets. 
%   Inputs:
%       data                data set (subjects x timepoints x regions of interest)
%       fdrV                fdr cutoff to be used
%       varargin    (optional, if a model is to be assumed)
%                   (both mmtype and lwflag have to be given as input)
%                   (use do_MM2 to get a mixture onto a predefined partial
%                       correation matrix)
%           usage:  do_MM(data,fdrV,mmtype,lwflag)
%           mmtype          mixture type to be fitted via EM
%                   'GGM'   Gauss-Gamma-Mixture
%                   'GIM'	Gauss-Inverse Gamma Mixture
%                   'LGM'	Laplace-Gamma Mixture
%                   'LIM'	Laplace-Inverse Gamma Mixture
%           lwFlag          1 if Ledoit-Wolf preprocessing to be used, else: 0
%       If no mmtype and lwFlag is given, optimal model is chosen based on Bayesian
%       Information Criterion (BIC) 
%
%   Outputs:
%   	res_MM              result mixture modeling 
%   	minBIC              BIC value of lowest result
%
% Author: Fabian Walocha (2017), fawalocha@gmail.com

%% Define important parameters

subjs = size(data,1);
ts = size(data,2);
roi = size(data,3);
connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
conns = sum(sum(connM));


%% Normalization of data

for IDX = 1:subjs
    tmp = data(IDX,:,:);
    data(IDX,:,:) = (tmp-mean(tmp(:)))./std(tmp(:));
end

%% Calculation inverse covariance matrix

[inv_covs1,inv_covs2] = deal(zeros([subjs,roi,roi]));

    
% Precision matrix to partial correlation matrix
for IDX1 = 1:subjs
    [sigma_hat,~] = covCor(squeeze(data(IDX1,:,:)));
    inv_covs1(IDX1,:,:) = inv(cov(squeeze(data(IDX1,:,:))));
    inv_covs2(IDX1,:,:) = inv(sigma_hat);
    for IDX2 = 1:roi
        for IDX3 = IDX2+1:roi
            inv_covs1(IDX1,IDX2,IDX3) = -inv_covs1(IDX1,IDX2,IDX3) / ...
                sqrt(inv_covs1(IDX1,IDX2,IDX2)*inv_covs1(IDX1,IDX3,IDX3));
            inv_covs1(IDX1,IDX3,IDX2) = -inv_covs1(IDX1,IDX3,IDX2) / ...
                sqrt(inv_covs1(IDX1,IDX2,IDX2)*inv_covs1(IDX1,IDX3,IDX3));
            inv_covs2(IDX1,IDX2,IDX3) = -inv_covs2(IDX1,IDX2,IDX3) / ...
                sqrt(inv_covs2(IDX1,IDX2,IDX2)*inv_covs2(IDX1,IDX3,IDX3));
            inv_covs2(IDX1,IDX3,IDX2) = -inv_covs2(IDX1,IDX3,IDX2) / ...
                sqrt(inv_covs2(IDX1,IDX2,IDX2)*inv_covs2(IDX1,IDX3,IDX3));
        end
    end
end

% If no model is demanded, a model is chosen via BIC
if nargin==2
    bics(1,:) = do_BIC(inv_covs1,'GGM');
    bics(2,:) = do_BIC(inv_covs2,'GGM');
    bics(3,:) = do_BIC(inv_covs1,'GIM');
    bics(4,:) = do_BIC(inv_covs2,'GIM');
    bics(5,:) = do_BIC(inv_covs1,'LGM');
    bics(6,:) = do_BIC(inv_covs2,'LGM');
    bics(7,:) = do_BIC(inv_covs1,'LIM');
    bics(8,:) = do_BIC(inv_covs2,'LIM');
    [minBIC,ind] = min(nanmean(bics,2));
    disp(minBIC);
    switch ind
        case 1
            inv_covs = inv_covs1;
            mmtype = 'GGM';
        case 2
            inv_covs = inv_covs2;
            mmtype = 'GGM';
        case 3
            inv_covs = inv_covs1;
            mmtype = 'GIM';
        case 4
            inv_covs = inv_covs2;
            mmtype = 'GIM';
        case 5
            inv_covs = inv_covs1;
            mmtype = 'LGM';
        case 6
            inv_covs = inv_covs2;
            mmtype = 'LGM';
        case 7
            inv_covs = inv_covs1;
            mmtype = 'LIM';
        case 8
            inv_covs = inv_covs2;
            mmtype = 'LIM';
    end
% if model criteria are given in varargin
else
    if lwFlag ==1
        inv_covs = inv_covs2;
    else
        inv_covs = inv_covs1;
    end
    mmtype = varargin{1};
end

res_MM = do_MM2(inv_covs,fdrV,mmtype);
 
end

