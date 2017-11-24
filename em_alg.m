function [params] = em_alg(data,mmtype,epsilon)
% EM_ALG Does mixture fit using EM algorithm
%   Input: 
%       data: partial correlation matrix for single subject
%       mmtype: 
%           'GGM' for Gauss-Gamma
%           'GIM' for Gauss-Inverse Gamma
%           'LGM' for Laplace-Gamma
%           'LIM' for Laplace-Inverse Gamma
%       epsilon: difference criterion over which the algorithm iterates
% Author: Fabian Walocha (2017), fawalocha@gmail.com

% InverseGamma or Gamma
if ismember('I',mmtype)
    params{1}.func = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
else
    params{1}.func = @gampdf;
end

% Laplace or Gauss
if ismember('L',mmtype)
    laplacepdf = @(x,mu,b) 1/(2*b) .* exp(-(abs(x-mu)./b));
    params{2}.func = laplacepdf;
else
    params{2}.func = @normpdf;
end

expectation = zeros([2,length(data)]);
chang = ones([2,1])*inf;

%% Initialization


%Initialize Gamma/Inverse Gamma
params{1}.p1 = 3;
params{1}.p2 = 1;
params{1}.alpha = 1/2;

% Initialize Gaussian/Laplace
params{2}.p1 = 0;
params{2}.p2 = 1;
params{2}.alpha = 1/2;


%% EM algorithm
while chang > epsilon
 
%% E step
locs = find(data<0);

expectation(1,:) = params{1}.alpha .* params{1}.func(data,params{1}.p1,params{1}.p2);
expectation(1,locs) = 0;
expectation(2,:) = params{2}.alpha  .* params{2}.func(data,params{2}.p1,params{2}.p2);        

% Normalize the expectation
expectation = expectation./repmat(sum(expectation,1),length(params),1);

%% M-step

for IDX = 1:2
    % Signal
    if IDX == 1 
        if strcmp(func2str(params{1}.func),'gampdf')
            % Gamma component   
            % E[X] = k * theta
            % Var[X] = k * theta²
            % theta = Var[X]/E[X]
            % k     = E[X]²/Var[X]
            mu = expectation(IDX,:)*data' / sum(expectation(IDX,:));
            var = expectation(IDX,:)*((data-mu).^2)' / sum(expectation(IDX,:));
            alpha = sum(expectation(IDX,:))/sum(expectation(:));
            if isnan(alpha)
                alpha = 10^(-1000);
            end
            chang = abs(params{IDX}.p1-mu^2/var) + abs(params{IDX}.p2-var/mu);
            params{IDX}.p2 = var/mu;
            params{IDX}.p1 = mu^2/var;
        else
            % Inverse Gamma component
            % E[X] = beta/(alpha-1)                     for alpha>1
            % Var[X] = beta^2 / ((alpha-1)^2(alpha-2))  for alpha>2
            % alpha = E[X]^2/Var[X]+2
            % beta  = E[X](E[X]^2/Var[X] +1)
            mu = expectation(IDX,:)*data' / sum(expectation(IDX,:));
            var = expectation(IDX,:)*((data-mu).^2)' / sum(expectation(IDX,:));
            alpha = sum(expectation(IDX,:))/sum(expectation(:));
            if isnan(alpha)
                alpha = 10^(-1000);
            end
            chang = abs(params{IDX}.p1-(mu^2/var+2)) + abs(params{IDX}.p2-(mu*(mu^2/var+2)));
            params{IDX}.p2 = mu*(mu^2/var+2);
            params{IDX}.p2(isnan(params{IDX}.p2)) = 0;
            params{IDX}.p1 = mu^2/var+2;
            params{IDX}.p1(isnan(params{IDX}.p1)) = 0;
        end
    % Noise
    elseif strcmp(func2str(params{IDX}.func),'normpdf')
        %Gauss component
        mu = expectation(IDX,:)*data' / sum(expectation(IDX,:));
        var = expectation(IDX,:)*((data-mu).^2)' / sum(expectation(IDX,:));  
        alpha = sum(expectation(IDX,:))/sum(expectation(:));
        if isnan(alpha)
            alpha = 10^(-1000);
        end
        chang = chang + abs(params{IDX}.p2-sqrt(var));
        % Gauss/Laplace mean parameter is not updated
        params{IDX}.p2 = sqrt(var);
    else
        %Laplace component
        mu = expectation(IDX,:)*data' / sum(expectation(IDX,:));
        var = expectation(IDX,:)*((data-mu).^2)' / sum(expectation(IDX,:));
        alpha = sum(expectation(IDX,:))/sum(expectation(:));
        if isnan(alpha)
            alpha = 10^(-1000);
        end
        chang = chang + abs(params{IDX}.p2-sqrt(var/2));
        % Gauss/Laplace mean parameter is not updated
        params{IDX}.p2 = sqrt(var/2);
        params{IDX}.p2(isnan(params{IDX}.p2)) = 0;
        
    end
    params{IDX}.alpha = alpha;
end    

end

end