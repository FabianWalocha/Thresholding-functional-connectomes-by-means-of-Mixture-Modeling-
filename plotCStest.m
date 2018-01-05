function plotCStest
% PLOTCSTEST Produces comparison of subject-wise vs cohort level methods
%   Applies both subject-wise and cohort-level permutation testing and
%   mixture modeling.
% Saves codes for: Figure 11 for the paper:
% Thresholding functional connectomes by means of Mixture Modeling (2017)  
% by Bielczyk, Walocha et al.
% Author: Fabian Walocha (2017), fawalocha@gmail.com

%% Prepare data

clc;
clearvars;
close all;

%% Set parameters

% real data or simulation data
dat = 'real';
mmtype = 'GGM';
lwFlag = 1;

% cutoff fdr value 
fdrV = 0.05;

for IDX11 = 1:2
    if IDX11 == 1 
        %% Real data (rest)

        cort_sys = 'visual';
        dataset = 'rest';

        pwdata = strcat('/data/HCP_',cort_sys,'_highpassfiltered_smoothed3mm_',dataset,'/');
        % Extract the previously identified 207 unrelated subjects
        tmt = csvread('data/SubjsParcUnr.txt');

        data = cell(length(tmt),1);
        for IDX = 1:length(tmt)
            data{IDX} = load(strcat(pwd,pwdata,cort_sys,num2str(tmt(IDX)),'.mat')); 
        end

        % Format data into SUBJSxHEMISPHERExTSxROI matrix
        data_mat = zeros([length(data),         2       ,length(data{1}.BOLD_left), size(data{1}.BOLD_left{1})]);
        for IDX1 = 1:length(data)
            for IDX2 = 1:length(data{1}.BOLD_left)
                for IDX3 = 1:size(data{1}.BOLD_left{1},2)
                   data_mat(IDX1,1,IDX2,:,IDX3) = ...
                                            data{IDX1}.BOLD_left{IDX2}(:,IDX3);
                   data_mat(IDX1,2,IDX2,:,IDX3) = ...
                                           data{IDX1}.BOLD_right{IDX2}(:,IDX3);
                end
            end
        end

        % Remove areas who don't survive parcellation and concatenate both
        % scanning directions
        data_mat_short = data_mat(:,:,:,:,[1:11 14:20]);

        % Concatenate both hemispheres
        d_short = cat(5,data_mat_short(:,1,:,:,:),data_mat_short(:,2,:,:,:));
        d_short = squeeze(d_short);

        % Reduce TS length to same length as Task data TS by concatenating both
        % scanning directions
        d_short = d_short(:,:,[6:205 1206:1405],:);
        data = squeeze(cat(3,d_short(:,1,:,:),d_short(:,2,:,:)));
    else
        %% Real data (task)

        cort_sys = 'visual';
        dataset = 'task';

        pwdata = strcat('/data/HCP_',cort_sys,'_highpassfiltered_smoothed3mm_',dataset,'/');
        % Extract the previously identified 207 unrelated subjects
        tmt = csvread('data/SubjsParcUnr.txt');

        data = cell(length(tmt),1);
        for IDX = 1:length(tmt)
            data{IDX} = load(strcat(pwd,pwdata,cort_sys,num2str(tmt(IDX)),'.mat')); 
        end

        % Format data into SUBJSxHEMISPHERExTSxROI matrix
        data_mat = zeros([length(data),         2       ,length(data{1}.BOLD_left), size(data{1}.BOLD_left{1})]);
         for IDX1 = 1:length(data)
             for IDX2 = 1:length(data{1}.BOLD_left)
                 for IDX3 = 1:size(data{1}.BOLD_left{1},2)
                    data_mat(IDX1,1,IDX2,:,IDX3) = ...
                                            data{IDX1}.BOLD_left{IDX2}(:,IDX3);
                    data_mat(IDX1,2,IDX2,:,IDX3) = ...
                                            data{IDX1}.BOLD_right{IDX2}(:,IDX3);
                 end
             end
         end

        % Remove areas who don't survive parcellation
        data_mat_short = data_mat(:,:,:,:,[1:11 14:20]);

        % Concatenate both hemispheres
        d_short = cat(5,data_mat_short(:,1,:,:,:),data_mat_short(:,2,:,:,:));
        d_short = squeeze(d_short);
        d_short_tmp = d_short(:,:,6:405,:);
        d_short = d_short(:,:,6:405,:);

        %concatenate both scanning directions
        d_short_tmp(:,1,:,:) = cat(3,d_short(:,1,[1:200],:),d_short(:,2,[1:200],:));
        d_short_tmp(:,2,:,:) = cat(3,d_short(:,1,[201:400],:),d_short(:,2,[201:400],:));
        d_short = d_short_tmp;
        data = squeeze(cat(3,d_short(:,1,:,:),d_short(:,2,:,:)));

    end

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

    %% Cohort Level Permutation Testing

    maxShuffle = 1000;
    timeseries = zeros([ts,roi]);
    res_permT2 = zeros([subjs,roi,roi]);
    for IDX10 = 1:subjs
        perm_covs = zeros([maxShuffle,roi,roi]);
        for IDX = 1:maxShuffle    
            for IDX2 = 1:roi
                timeseries(:,IDX2) = data(IDX10,randperm(ts),IDX2);
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
        for IDX2 = 1:roi
            for IDX3 = 1:roi
                res_permT2(IDX10,IDX2,IDX3) = ...
                    sum(perm_covs(:,IDX2,IDX3)>par_corrs1(IDX10,IDX2,IDX3))/maxShuffle;       
            end
        end
    end

    %% Single subject Permutation testing

    % Permutation testing
    maxShuffle = 1000;
    p=0.05;
    perm_covs2 = zeros([maxShuffle,roi,roi]);
    % Get 1000 partial correlation matrices where time series is shuffled
    % for each set of ROIs
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
        perm_covs2(IDX,:,:) = precis;
    end

    res_permT1 = zeros([subjs,roi,roi]);
    for IDX1 = 1:subjs
        for IDX2 = 1:roi
            for IDX3 = 1:roi
                res_permT1(IDX1,IDX2,IDX3) = ...
                    (sum(perm_covs2(:,IDX2,IDX3)>par_corrs1(IDX1,IDX2,IDX3))/maxShuffle)<p;       
            end
        end
    end

    %% Single Subject 

    thresh = zeros(subjs,1);
    for IDX = 1:subjs

        tmp = squeeze(par_corrs(IDX,:,:));
        tmp = tmp(connM);
        datavec = tmp;

        dataMM = (datavec./std(datavec))';

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

    res_MM1 = double(par_corrs>=thresh);

    %% Cohort Level MM


    datavec = zeros([subjs,conns]);
    for IDX = 1:subjs
        tmp = squeeze(par_corrs(IDX,:,:));
        tmp = tmp(connM);
        datavec(IDX,:) = tmp;
    end

    datavec = squeeze(reshape(datavec,[subjs*conns,1]));
    dataMM2 = (datavec./std(datavec))';

    params = em_alg(dataMM2,mmtype,0.001);

    [datS,~] = sort(dataMM2);
    
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
        thresh = +inf;
    else
        thresh = datS(x(1));
    end
    thresh = thresh * std(datavec);

    res_MM2 = double(par_corrs>=thresh);

    if IDX11 == 1
        res_permT1rest = res_permT1;
        res_permT2rest = res_permT2;
        res_MM1rest = res_MM1;
        res_MM2rest = res_MM2;
    else
        res_permT1task = res_permT1;
        res_permT2task = res_permT2;
        res_MM1task = res_MM1;
        res_MM2task = res_MM2;
    end
end

nams = {'permT1','permT2','MM1','MM2'};

for IDX = 1:4
    switch IDX
        case 1
            resultR = res_permT1rest;
            resultT = res_permT1task;
        case 2
            resultR = res_permT2rest;
            resultT = res_permT2task;
        case 3
            resultR = res_MM1rest;
            resultT = res_MM1task;
        case 4 
            resultR = res_MM2rest;
            resultT = res_MM2task;
    end
    
    resultD = squeeze(sum(resultT-resultR,1));
    mwu2 = mannwhitu(resultR,resultT);
    resultD3 = resultD.*mwu2;
   
    switch IDX
        case 1
            res_permT1 = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_permT1');

            res_permT1 = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_permT1');
            
        case 2
            res_permT2 = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_permT2');

            res_permT2 = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_permT2');
        case 3
            res_MM1 = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_MM1');

            res_MM1 = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_MM1');
        case 4
            res_MM2 = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_MM2');

            res_MM2 = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_MM2');
    end


end

res_permT1 = squeeze(sum(res_permT1rest,1));
res_permT2 = squeeze(sum(res_permT2rest,1));
res_MM1 = squeeze(sum(res_MM1rest,1));
res_MM2 = squeeze(sum(res_MM2rest,1));
save('results/res_permT1rest','res_permT1');
save('results/res_permT2rest','res_permT2');
save('results/res_MM1rest','res_MM1');
save('results/res_MM2rest','res_MM2');

res_permT1 = squeeze(sum(res_permT1task,1));
res_permT2 = squeeze(sum(res_permT2task,1));
res_MM1 = squeeze(sum(res_MM1task,1));
res_MM2 = squeeze(sum(res_MM2task,1));
save('results/res_permT1task','res_permT1');
save('results/res_permT2task','res_permT2');
save('results/res_MM1task','res_MM1');
save('results/res_MM2task','res_MM2');    

end

function [ mwu ] = mannwhitu( matR,matT )
%MANNWHITU performs Mann-Whitney-U test with Bonferroni correction
%connection-wise between rest and task
%   Input:
%           matR = connectivity matrices in rest 
%           matT = connectivity matrices in task
%
%   Output:
%           mwu = Significant connections via ranksum


roi = size(matR,2); 
connM = triu(ones(roi))-diag(diag(ones(roi)));
conns = sum(sum(connM));
subjs = size(matR,1);

mwu = zeros([roi,roi]);
for IDX1 = 1:roi
    for IDX2 = IDX1+1:roi
        mwu(IDX2,IDX1) = ranksum(matR(:,IDX1,IDX2),matT(:,IDX1,IDX2))<(0.01/conns);
        mwu(IDX1,IDX2) = ranksum(matR(:,IDX1,IDX2),matT(:,IDX1,IDX2))<(0.01/conns);
    end
end

end