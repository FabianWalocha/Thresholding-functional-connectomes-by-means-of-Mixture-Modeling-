function [ res_MM,minBIC ] = do_MM2( par_corrs,fdrV,varargin )
% DO_MM2 Employs mixture modeling based approach on a set of partial
% correlation matrices
%   Inputs:
%           par_corrs       set of partial correlation matrices
%           fdrV            cutoff FDR value to be employed
%           varargin    (optional, if a model is to be assumed)
%               mmtype  mixture type to be fitted via EM
%               	'GGM'   Gauss-Gamma-Mixture
%                   'GIM'	Gauss-Inverse Gamma Mixture
%                   'LGM'	Laplace-Gamma Mixture
%                   'LIM'	Laplace-Inverse Gamma Mixture
%               If no mmtype and lwFlag is given, optimal model is chosen based on 
%               Bayesian Information Criterion (BIC) 
%   Outputs:
%           res_MM          result mixture modeling 
%           minBIC          BIC value of lowest result
% Author: Fabian Walocha (2017), fawalocha@gmail.com

subjs = size(par_corrs,1);
roi = size(par_corrs,2);
connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
conns = sum(sum(connM));

if nargin==2
    bics(1,:) = do_BIC(par_corrs1,'GGM');
    bics(2,:) = do_BIC(par_corrs1,'GIM');
    bics(3,:) = do_BIC(par_corrs1,'LGM');
    bics(4,:) = do_BIC(par_corrs1,'LIM');
    [minBIC,ind] = min(nanmean(bics,2));
    switch ind
        case 1
            mmtype = 'GGM';
        case 2
            mmtype = 'GIM';
        case 3
            par_corrs = par_corrs2;
            mmtype = 'LGM';
        case 4
            par_corrs = par_corrs1;
            mmtype = 'LIM';
    end
else
% if model criteria are given in varargin
    mmtype = varargin{1};
end
%----------------------- Part 2 ------------------------------------
thresh = zeros(subjs,1);
for IDX = 1:subjs

    tmp = squeeze(par_corrs(IDX,:,:));
    tmp = tmp(connM);
    datavec = tmp;

    dataMM = (datavec./std(datavec))';
    
    % Estimating FDR by means of EM
    params = em_alg(dataMM,mmtype,0.001);

    [datS,~] = sort(dataMM);
    
    
    if strcmp(mmtype,'LIM') || strcmp(mmtype,'LGM')
        laplacecdf = @(x,mu,b) (x<mu) .* 0.5.*exp((x-mu)./b) + (x>=mu) .* (1-0.5.*exp(-(x-mu)./b));
        fpr = params{2}.alpha .* (1-(laplacecdf(datS,params{2}.p1,params{2}.p2)));
    else 
        fpr = params{2}.alpha .* (1-(normcdf(datS,params{2}.p1,params{2}.p2)));
    end

    ppr = 1-((1:length(datS))./length(datS));

    fdr = fpr./ppr;
    
    % return all indices for values bigger than thresh
    x = find(fdr<fdrV);

    % If FDR cannot be determined, set threshold to +infinty
    if isempty(x)
        thresh(IDX) = +inf;
    else
        thresh(IDX) = datS(x(1));
    end
    thresh(IDX) = thresh(IDX) * std(datavec);
end

% Thresholding according to  for each partial correlation matrix
thresh = repmat(thresh,[1,roi,roi]);
    
res_MM = par_corrs>=thresh;

end

