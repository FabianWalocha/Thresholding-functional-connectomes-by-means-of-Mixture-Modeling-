function [ res_MM, res_permT, res_EP, res_ledW, res_N1, res_N2, iccs, avgs, icc2 ] = plotMain1( data,fdrV,lwFlag,rOrT,mmtype,varargin)
%PLOTMAIN1 Applies various methods of sparsifying connectomes onto data
%   Inputs:
%       data        data of type subjsxtimepointsxrois
%       mmtype      
%                   'GGM' for Gauss-Gamma
%                   'GIM' for Gauss-InverseGamma
%                   'LGM' for Laplace-Gamma
%                   'LIM' for Laplace-InverseGamma 
%                   
%       fdrV        ciritcal false discovery value (threshold)
%       lwFlag      Use Ledoit-Wolf regularization?
%               
%       rOrT        rest or task? (used only for saving data and plotting)
%                   0 = sim, 1 = rest, 2 = task
%       varargin
%           valMat  validation matrix if data is simulated data
%   Outputs:
%       res_MM      result mixture modeling 
%       res_permT   result permutation testing 
%       res_EP      result empirical precision, thresholded at 0
%       res_ledW    result ledoit-wolf, thresholded at 0
%       res_N1      result thresholding upper 5%
%       res_N2      result thresholding upper 10%
%       iccs        Index Overlap
%       avgs        averages
%       icc2        cohort-level ICC 
%
% Produces: Figure 1A, figure 1B, Figure 3, table 1, table 2, table 3,
% table 4, saves results for figure 9, figure 10 for the paper:
% Thresholding functional connectomes by means of Mixture Modeling (2017)  
% by Bielczyk, Walocha et al.
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

[par_corrs1,par_corrs2] = deal(zeros([subjs,roi,roi]));

    
% Precision values to partial correlations
for IDX1 = 1:subjs
    [sigma_hat,~] = covCor(squeeze(data(IDX1,:,:)));
    par_corrs1(IDX1,:,:) = inv(cov(squeeze(data(IDX1,:,:))));
    par_corrs2(IDX1,:,:) = inv(sigma_hat);
    for IDX2 = 1:roi
        for IDX3 = IDX2+1:roi
            par_corrs1(IDX1,IDX2,IDX3) = -par_corrs1(IDX1,IDX2,IDX3) / ...
                sqrt(par_corrs1(IDX1,IDX2,IDX2)*par_corrs1(IDX1,IDX3,IDX3));
            par_corrs1(IDX1,IDX3,IDX2) = -par_corrs1(IDX1,IDX3,IDX2) / ...
                sqrt(par_corrs1(IDX1,IDX2,IDX2)*par_corrs1(IDX1,IDX3,IDX3));
            par_corrs2(IDX1,IDX2,IDX3) = -par_corrs2(IDX1,IDX2,IDX3) / ...
                sqrt(par_corrs2(IDX1,IDX2,IDX2)*par_corrs2(IDX1,IDX3,IDX3));
            par_corrs2(IDX1,IDX3,IDX2) = -par_corrs2(IDX1,IDX3,IDX2) / ...
                sqrt(par_corrs2(IDX1,IDX2,IDX2)*par_corrs2(IDX1,IDX3,IDX3));
        end
    end
end

if lwFlag ==1
    par_corrs = par_corrs2;
else
    par_corrs = parr_corrs1;
end


subjs = size(par_corrs,1);
roi = size(par_corrs,2);
connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
conns = sum(sum(connM));


thresh = zeros(subjs,1);
for IDX = 1:subjs

    tmp = squeeze(par_corrs(IDX,:,:));
    tmp = tmp(connM);
    datavec = tmp;

    dataMM = (datavec./std(datavec))';

    params = em_alg(dataMM,mmtype,0.001);

    [datS,~] = sort(dataMM);
    if IDX == 53
        fig1 = figure('visible','off');
        dataP = dataMM(dataMM>-4);
        dataP = dataP(dataP<6);
        [f,x] = hist(dataP,50);
        bar(x,f/trapz(x,f));
        hold on

        range=-4:.01:6;

        plot(range,params{1}.alpha.*params{1}.func(range,params{1}.p1,params{1}.p2),'g','LineWidth',2);
        plot(range,params{2}.alpha.*params{2}.func(range,params{2}.p1,params{2}.p2),'r','LineWidth',2);

        ylabel('probability density','Fontsize',12)
        ylim([0,1])

        legend('Partial correlation values','Est. signal distribution','Est. noise distribution','Location','nw')
        if rOrT == 1
            saveTightFigure(fig1,strcat('plots/modelFit_',mmtype,'_rest'));
        elseif rOrT == 2
            saveTightFigure(fig1,strcat('plots/modelFit_',mmtype,'_task'));                
        end
    end

    if strcmp(mmtype,'LIM') || strcmp(mmtype,'LGM')
        laplacecdf = @(x,mu,b) (x<mu) .* 0.5.*exp((x-mu)./b) + (x>=mu) .* (1-0.5.*exp(-(x-mu)./b));
        fpr = params{2}.alpha .* (1-(laplacecdf(datS,params{2}.p1,params{2}.p2)));
    else 
        fpr = params{2}.alpha .* (1-(normcdf(datS,params{2}.p1,params{2}.p2)));

    end

    ppr = 1-((1:length(datS))./length(datS));

    fdr = fpr./ppr;

    % Automatic thresolding at the cross-section

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

% One threshold for each partial correlation matrix
thresh = repmat(thresh,[1,roi,roi]);
    
res_MM = par_corrs>=thresh;

%% Comparison to other methods

% Empirical precision, thresholded at 0
res_EP = par_corrs1>0;


% Ledoit-Wolf regularized, thresholded at 0
res_ledW = par_corrs2>0;

% Permutation testing
maxShuffle = 1000;
p=0.05;
perm_covs = zeros([maxShuffle,roi,roi]);
for IDX = 1:maxShuffle    
    rp = randperm(subjs);
    rp = rp(1:roi);
    timeseries = zeros([ts,roi]);
    for IDX2 = 1:roi
        timeseries(:,IDX2) = data(rp(IDX2),:,IDX2);
    end
    precis = inv(cov(timeseries));
    for IDX1 = 1:roi
        for IDX2 = IDX1+1:roi
            precis(IDX1,IDX2) = -precis(IDX1,IDX2)/sqrt(precis(IDX1,IDX1)*precis(IDX2,IDX2));
            precis(IDX2,IDX1) = -precis(IDX2,IDX1)/sqrt(precis(IDX1,IDX1)*precis(IDX2,IDX2));
        end
    end
    perm_covs(IDX,:,:) = precis;
end

res_permT = zeros([subjs,roi,roi]);
for IDX1 = 1:subjs
    for IDX2 = 1:roi
        for IDX3 = 1:roi
            res_permT(IDX1,IDX2,IDX3) = ...
                (sum(perm_covs(:,IDX2,IDX3)>par_corrs1(IDX1,IDX2,IDX3))/maxShuffle)<p;       
        end
    end
end

[res_N1,res_N2] = deal(par_corrs1);
for IDX = 1:subjs
    tmp = squeeze(par_corrs1(IDX,:,:));
    tmp = sort(tmp(connM),'descend');
    res_N1(IDX,:,:) = par_corrs1(IDX,:,:) > tmp(floor(conns*0.05));
    res_N2(IDX,:,:) = par_corrs1(IDX,:,:) > tmp(floor(conns*0.1));
end

%% Comparing the performance (if simulation)

if nargin>5
    net = varargin{1};  
    net = net ~= 0;

    [r1,r2,r3,r4, r5, r6, tpr1,tpr2,tpr3,tpr4,tpr5,tpr6,fpr1,fpr2,fpr3,fpr4,fpr5,fpr6] = deal(zeros(subjs,1));

    for IDX = 1:subjs
        tmp1 = squeeze(res_EP(IDX,:,:));
        tmp2 = squeeze(res_ledW(IDX,:,:));
        tmp3 = squeeze(res_permT(IDX,:,:));
        tmp4 = squeeze(res_MM(IDX,:,:));
        tmp5 = squeeze(res_N1(IDX,:,:));
        tmp6 = squeeze(res_N2(IDX,:,:));
        tmpnet = squeeze(net(IDX,:,:));

        tmp1 = tmp1(connM);
        tmp2 = tmp2(connM);
        tmp3 = tmp3(connM);
        tmp4 = tmp4(connM);
        tmp5 = tmp5(connM);
        tmp6 = tmp6(connM);
        tmpnet = tmpnet(connM);

        tpr1(IDX) = sum(tmp1.*tmpnet)/sum(tmpnet~=0);
        tpr2(IDX) = sum(tmp2.*tmpnet)/sum(tmpnet~=0);
        tpr3(IDX) = sum(tmp3.*tmpnet)/sum(tmpnet~=0);
        tpr4(IDX) = sum(tmp4.*tmpnet)/sum(tmpnet~=0);
        tpr5(IDX) = sum(tmp5.*tmpnet)/sum(tmpnet~=0);
        tpr6(IDX) = sum(tmp6.*tmpnet)/sum(tmpnet~=0);

        fpr1(IDX) = sum(tmp1.*(tmpnet==0))/sum(tmpnet==0);
        fpr2(IDX) = sum(tmp2.*(tmpnet==0))/sum(tmpnet==0);
        fpr3(IDX) = sum(tmp3.*(tmpnet==0))/sum(tmpnet==0);
        fpr4(IDX) = sum(tmp4.*(tmpnet==0))/sum(tmpnet==0);
        fpr5(IDX) = sum(tmp5.*(tmpnet==0))/sum(tmpnet==0);
        fpr6(IDX) = sum(tmp6.*(tmpnet==0))/sum(tmpnet==0);

        r1(IDX) = sum(tmp1==tmpnet)/conns;
        r2(IDX) = sum(tmp2==tmpnet)/conns;
        r3(IDX) = sum(tmp3==tmpnet)/conns;
        r4(IDX) = sum(tmp4==tmpnet)/conns;
        r5(IDX) = sum(tmp5==tmpnet)/conns;
        r6(IDX) = sum(tmp6==tmpnet)/conns;
    end


    perf(1,1) = mean(r1);
    perf(2,1) = mean(r2);
    perf(3,1) = mean(r3);
    perf(4,1) = mean(r4);
    perf(5,1) = mean(r5);
    perf(6,1) = mean(r6);
    
    perf(1,2) = std(r1);
    perf(2,2) = std(r2);
    perf(3,2) = std(r3);
    perf(4,2) = std(r4);
    perf(5,2) = std(r5);
    perf(6,2) = std(r6);    
        
    tp(1,1) = mean(tpr1);
    tp(2,1) = mean(tpr2);
    tp(3,1) = mean(tpr3);
    tp(4,1) = mean(tpr4);
    tp(5,1) = mean(tpr5);
    tp(6,1) = mean(tpr6);
    
    tp(1,2) = std(tpr1);
    tp(2,2) = std(tpr2);
    tp(3,2) = std(tpr3);
    tp(4,2) = std(tpr4);
    tp(5,2) = std(tpr5);
    tp(6,2) = std(tpr6);
    
    fp(1,1) = mean(fpr1);
    fp(2,1) = mean(fpr2);
    fp(3,1) = mean(fpr3);
    fp(4,1) = mean(fpr4);
    fp(5,1) = mean(fpr5);
    fp(6,1) = mean(fpr6);
    
    fp(1,2) = std(fpr1);
    fp(2,2) = std(fpr2);
    fp(3,2) = std(fpr3);
    fp(4,2) = std(fpr4);
    fp(5,2) = std(fpr5);
    fp(6,2) = std(fpr6);    
    
    model_series = [perf(:,1),tp(:,1),fp(:,1)];
    model_error = [perf(:,2),tp(:,2),fp(:,2)];
    ax = axes;
    fig2 = figure('visible','off');
    bar(model_series,'BarWidth',1.5);
     xticks(ax,[1 2 3 4 5]);
    hold on;
    ngroups = size(model_series, 1);
    nbars = size(model_series, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end

    legend('Mean performance','TPR','FPR','Location','southwest');
    names = {'Empirical precision','Ledoit-Wolf','Permutation Testing','Mixture modeling','Prop. thresholding, 5%','Prop. thresholding, 10%'};
    set(gca,'XTickLabel',names,'fontsize',4)

    saveTightFigure(fig2,'plots/perfSim');
    
else
    
    %% Means of each method
    avgs(1) = mean(sum(res_EP(:,connM),2))/conns;
    avgs(2) = mean(sum(res_ledW(:,connM),2))/conns;
    avgs(3) = mean(sum(res_permT(:,connM),2))/conns;
    avgs(4) = mean(sum(res_MM(:,connM),2))/conns;
    avgs(5) = mean(sum(res_N1(:,connM),2))/conns;
    avgs(6) = mean(sum(res_N2(:,connM),2))/conns;
    avgs(7) = std(sum(res_EP(:,connM),2))/conns;
    avgs(8) = std(sum(res_ledW(:,connM),2))/conns;
    avgs(9) = std(sum(res_permT(:,connM),2))/conns;
    avgs(10) = std(sum(res_MM(:,connM),2))/conns;
    avgs(11) = std(sum(res_N1(:,connM),2))/conns;
    avgs(12) = std(sum(res_N2(:,connM),2))/conns;
    
  
    
    %% Figure 4: Check test-retest reliability (on real data)
    
    dat1 = data(:,(1:floor(end/2)),:);
    dat2 = data(:,((floor(end/2)+1):(floor(end/2)*2)),:);
    
    ts = ts/2;
    
    icc2 = zeros([6,1]);
    
    for IDX4 = 1:2
        
        if IDX4 == 1
            data = dat1;
        else
            data = dat2;
        end
        %% Calculation inverse covariance matrix

        [par_corrs1,par_corrs2] = deal(zeros([subjs,roi,roi]));
 

        for IDX1 = 1:subjs
            [sigma_hat,~] = covCor(squeeze(data(IDX1,:,:)));
            par_corrs1(IDX1,:,:) = inv(cov(squeeze(data(IDX1,:,:))));
            par_corrs2(IDX1,:,:) = inv(sigma_hat);
            for IDX2 = 1:size(par_corrs1,2)
                for IDX3 = IDX2+1:size(par_corrs1,2)
                    par_corrs1(IDX1,IDX2,IDX3) = -par_corrs1(IDX1,IDX2,IDX3) / ...
                        sqrt(par_corrs1(IDX1,IDX2,IDX2)*par_corrs1(IDX1,IDX3,IDX3));
                    par_corrs1(IDX1,IDX3,IDX2) = -par_corrs1(IDX1,IDX3,IDX2) / ...
                        sqrt(par_corrs1(IDX1,IDX2,IDX2)*par_corrs1(IDX1,IDX3,IDX3));
                    par_corrs2(IDX1,IDX2,IDX3) = -par_corrs2(IDX1,IDX2,IDX3) / ...
                        sqrt(par_corrs2(IDX1,IDX2,IDX2)*par_corrs2(IDX1,IDX3,IDX3));
                    par_corrs2(IDX1,IDX3,IDX2) = -par_corrs2(IDX1,IDX3,IDX2) / ...
                        sqrt(par_corrs2(IDX1,IDX2,IDX2)*par_corrs2(IDX1,IDX3,IDX3));
                end
            end
        end
            
        if lwFlag 
            res_MMicc = do_MM2(par_corrs2,fdrV,mmtype);
        else
            res_MMicc = do_MM2(par_corrs1,fdrV,mmtype);
        end
        
        [res_N1icc,res_N2icc] = deal(par_corrs1);
        for IDX = 1:subjs
            tmp = squeeze(par_corrs1(IDX,:,:));
            tmp = sort(tmp(connM),'descend');
            res_N1icc(IDX,:,:) = par_corrs1(IDX,:,:) > tmp(floor(conns*0.05));
            res_N2icc(IDX,:,:) = par_corrs1(IDX,:,:) > tmp(floor(conns*0.1));
        end

        
        %% Comparison to other methods

        % Empirical precision, thresholded at 0
        res_EPicc = par_corrs1>0;

        % Ledoit-Wolf regularized, thresholded at 0
        res_ledWicc = par_corrs2>0;

        % Permutation testing
        maxShuffle = 1000;
        p=0.05;
        perm_covs = zeros([maxShuffle,roi,roi]);
        for IDX = 1:maxShuffle    
            rp = randperm(subjs);
            rp = rp(1:roi);
            timeseries = zeros([ts,roi]);
            for IDX2 = 1:roi
                timeseries(:,IDX2) = data(rp(IDX2),:,IDX2);
            end
            precis = inv(cov(timeseries));
            for IDX1 = 1:size(precis,1)
                for IDX2 = IDX1+1:size(precis,1)
                    precis(IDX1,IDX2) = -precis(IDX1,IDX2)/sqrt(precis(IDX1,IDX1)*precis(IDX2,IDX2));
                    precis(IDX2,IDX1) = -precis(IDX2,IDX1)/sqrt(precis(IDX1,IDX1)*precis(IDX2,IDX2));
                end
            end
            perm_covs(IDX,:,:) = precis;
        end

        res_permTicc = zeros([subjs,roi,roi]);
        for IDX1 = 1:subjs
            for IDX2 = 1:roi
                for IDX3 = 1:roi
                    res_permTicc(IDX1,IDX2,IDX3) = ...
                        sum(perm_covs(:,IDX2,IDX3)>par_corrs1(IDX1,IDX2,IDX3))/maxShuffle<p;       
                end
            end
        end
        

        if IDX4 == 1
            res{1,1} = res_EPicc;
            res{1,2} = res_ledWicc;
            res{1,3} = res_permTicc;
            res{1,4} = res_MMicc;
            res{1,5} = res_N1icc;
            res{1,6} = res_N2icc;
        else 
            res{2,1} = res_EPicc;
            res{2,2} = res_ledWicc;
            res{2,3} = res_permTicc;
            res{2,4} = res_MMicc;
            res{2,5} = res_N1icc;
            res{2,6} = res_N2icc;
        end
        if(IDX4 == 1)
            for IDX11 = 1:6
                icc2(IDX11) = ICC_new(3,'single',res{1,IDX11}(:,connM)');
            end
        end
    end

    res_ICC = cell(size(res,2),subjs);
    % Calculate ICC score subjectwise
    [iccs] = deal(zeros(size(res_ICC,1),subjs));
    for IDX1 = 1:size(res_ICC,1)
        res_ICC{IDX1} = [];
        for IDX2 = 1:subjs 
            tmp1 = squeeze(res{1,IDX1}(IDX2,:,:));
            tmp1 = tmp1(connM);
            tmp2 = squeeze(res{2,IDX1}(IDX2,:,:));
            tmp2 = tmp2(connM);
            res_ICC{IDX1,IDX2} = [tmp1(:),tmp2(:)];
            iccs(IDX1,IDX2) = sum(tmp1==tmp2)/conns;
        end
        
    end 
    
   
    
    % Save result matrices for use in circular graphs
    res_EP3 = squeeze(sum(res_EP,1));
    res_ledW3 = squeeze(sum(res_ledW,1));
    res_permT3 = squeeze(sum(res_permT,1));
    res_MM3 = squeeze(sum(res_MM,1));
    res_N13 = squeeze(sum(res_N1,1));
    res_N23 = squeeze(sum(res_N2,1));
    
    if rOrT == 1
        save('results/res_EPrest','res_EP3');
        save('results/res_ledWrest','res_ledW3');
        save('results/res_permTrest','res_permT3');
        save('results/res_MMrest','res_MM3');
        save('results/res_N1rest','res_N13');
        save('results/res_N2rest','res_N23');
    else
        save('results/res_EPtask','res_EP3');
        save('results/res_ledWtask','res_ledW3');
        save('results/res_permTtask','res_permT3');
        save('results/res_MMtask','res_MM3');
        save('results/res_N1task','res_N13');
        save('results/res_N2task','res_N23');
    end
end

