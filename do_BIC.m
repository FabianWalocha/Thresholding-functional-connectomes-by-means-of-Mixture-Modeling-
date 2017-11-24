function [ bics ] = do_BIC( inv_covs,mmtype )
%DO_MIXMOD2 Calculates the Bayesian Information Criterion for a given set
%of partial correlation matrices and a mixture type
%   Inputs:
%       inv_covs    set of given partial correlation matrices in the 
%                   form of (subjects x regions of interest x regions 
%                   of interest)
%       mmtype  mixture type to be fitted via EM
%       	'GGM'	Gauss-Gamma-Mixture
%           'GIM'	Gauss-Inverse Gamma Mixture
%           'LGM'	Laplace-Gamma Mixture
%           'LIM'	Laplace-Inverse Gamma Mixture
%   Outputs:
%       bics        Bayesian Information Criterion for the selected model
% Author: Fabian Walocha (2017), fawalocha@gmail.com

subjs = size(inv_covs,1);
roi = size(inv_covs,2);
connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
bics = zeros(subjs,1);

for IDX = 1:subjs

    tmp = squeeze(inv_covs(IDX,:,:));
    tmp = tmp(connM);
    datavec = tmp;

    dataMM = (datavec./std(datavec))';

    params = em_alg(dataMM,mmtype,0.001);

    locs = find(dataMM<0);
    a = params{1}.alpha .*params{1}.func(dataMM,params{1}.p1,params{1}.p2);
    a(locs) = 0;
    b = params{2}.alpha .*params{2}.func(dataMM,params{2}.p1,params{2}.p2);
    % Number of parameters is 4 for all used mixtures
    numPar = 4;
    
    bics(IDX) = numPar * log(length(dataMM))+(-2)*sum(log(a+b)); 
end
end