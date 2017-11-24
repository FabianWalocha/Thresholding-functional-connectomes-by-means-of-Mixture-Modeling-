function plotVarN
% PLOTVARN Tests the results for stability over decreasing number of observed ROIs
% Produces: Figure 2
% Author: Fabian Walocha (2017), fawalocha@gmail.com

clc;
clearvars;
close all;

%% Set parameters

% real data or simulation data
% dat = 'real';
dat = 'sim';
mmtype = 'GGM';
lwFlag = 1;

% cutoff fdr value 
fdrV = 0.05;

%% Simulation data

% Read data
data_mat = importdata(strcat(pwd,'/data/Smith_orig/sim4.mat'));

subjs = data_mat.Nsubjects;
roi = data_mat.Nnodes;
ts = size(data_mat.ts,1)/subjs;
net = data_mat.net;
netF = net ~= 0;
connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
conns = sum(sum(connM));


% Format data into SUBJSxTSxROI matrix
data = zeros([subjs,ts,roi]);
for IDX1 = 1:subjs
    data(IDX1,:,:) = real(data_mat.ts((1+ts*(IDX1-1)):(ts*IDX1),:));
end

%% Define important parameters

subjs = size(data,1);
ts = size(data,2);
roi = size(data,3);

%% Normalization of data

for IDX = 1:subjs
    tmp = data(IDX,:,:);
    data(IDX,:,:) = (tmp-mean(tmp(:)))./std(tmp(:));
end

%% Calculation inverse covariance matrix

[par_corrs1F,par_corrs2F] = deal(zeros([subjs,roi,roi]));


% Precision values to partial correlations
for IDX1 = 1:subjs
    [sigma_hat,~] = covCor(squeeze(data(IDX1,:,:)));
    par_corrs1F(IDX1,:,:) = inv(cov(squeeze(data(IDX1,:,:))));
    par_corrs2F(IDX1,:,:) = inv(sigma_hat);
    for IDX2 = 1:roi
        for IDX3 = IDX2+1:roi
            par_corrs1F(IDX1,IDX2,IDX3) = -par_corrs1F(IDX1,IDX2,IDX3) / ...
                sqrt(par_corrs1F(IDX1,IDX2,IDX2)*par_corrs1F(IDX1,IDX3,IDX3));
            par_corrs1F(IDX1,IDX3,IDX2) = -par_corrs1F(IDX1,IDX3,IDX2) / ...
                sqrt(par_corrs1F(IDX1,IDX2,IDX2)*par_corrs1F(IDX1,IDX3,IDX3));
            par_corrs2F(IDX1,IDX2,IDX3) = -par_corrs2F(IDX1,IDX2,IDX3) / ...
                sqrt(par_corrs2F(IDX1,IDX2,IDX2)*par_corrs2F(IDX1,IDX3,IDX3));
            par_corrs2F(IDX1,IDX3,IDX2) = -par_corrs2F(IDX1,IDX3,IDX2) / ...
                sqrt(par_corrs2F(IDX1,IDX2,IDX2)*par_corrs2F(IDX1,IDX3,IDX3));
        end
    end
end

par_corrs1 = par_corrs1F;
par_corrs2 = par_corrs2F;
net = netF;
revNs = zeros([45,1]);
[fp,tp,perf] = deal([6,45]);

for IDX10 = 1:45

    roi = roi-1;
    connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
    conns = sum(sum(connM));
    
    revN = randi(roi,1,1);
    revNs(IDX10) = revN;
    
    par_corrs1 = par_corrs1(:,1:end~=revN,1:end~=revN);
    par_corrs2 = par_corrs2(:,1:end~=revN,1:end~=revN);
    net = net(:,1:end~=revN,1:end~=revN);
    
    data = data(:,:,1:end~=revN); 
    

    if lwFlag ==1
        res_MM = do_MM2(par_corrs2,fdrV,mmtype);
    else
        res_MM = do_MM2(par_corrs1,fdrV,mmtype);
    end


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
    
    % Proportional thresholding (5%/10%)
    if conns*0.1<1
        res_N1 = par_corrs1>=inf;
        res_N2 = par_corrs1>=inf;
    else
        if conns*0.05 <1
            res_N1 = par_corrs1>=inf;
            res_N2 = par_corrs1;
            for IDX = 1:subjs
                tmp = squeeze(par_corrs1(IDX,:,:));
                tmp = sort(tmp(connM),'descend');
                res_N2(IDX,:,:) = par_corrs1(IDX,:,:) >= tmp(floor(conns*0.1));
            end
        else
            [res_N1,res_N2] = deal(par_corrs1);
            for IDX = 1:subjs
                tmp = squeeze(par_corrs1(IDX,:,:));
                tmp = sort(tmp(connM),'descend');

                res_N1(IDX,:,:) = par_corrs1(IDX,:,:) >= tmp(floor(conns*0.05));
                res_N2(IDX,:,:) = par_corrs1(IDX,:,:) >= tmp(floor(conns*0.1));
            end
        end
    end

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


    perf(1,IDX10) = mean(r1);
    perf(2,IDX10) = mean(r2);
    perf(3,IDX10) = mean(r3);
    perf(4,IDX10) = mean(r4);
    perf(5,IDX10) = mean(r5);
    perf(6,IDX10) = mean(r6);
    
    
    tp(1,IDX10) = nanmean(tpr1);
    tp(2,IDX10) = nanmean(tpr2);
    tp(3,IDX10) = nanmean(tpr3);
    tp(4,IDX10) = nanmean(tpr4); 
    tp(5,IDX10) = nanmean(tpr5);
    tp(6,IDX10) = nanmean(tpr6);
    
    fp(1,IDX10) = mean(fpr1);
    fp(2,IDX10) = mean(fpr2);
    fp(3,IDX10) = mean(fpr3);
    fp(4,IDX10) = mean(fpr4);
    fp(5,IDX10) = mean(fpr5);
    fp(6,IDX10) = mean(fpr6);

end

f1 =figure('visible','off');
plot(49:-1:5,perf,'--o','linewidth',1.5);
set(gca,'xdir','reverse')
legend({'Empirical precision','Ledoit-Wolf','Permutation Testing','Mixture modeling',...
    'Proportional thresholding, top 5%','Proportional thresholding, top 10%'},'location','southwest','FontSize',8)
xlabel('Number of nodes in network','FontSize',12)
ylabel('Mean performance','FontSize',12)

f2 = figure('visible','off');
plot(49:-1:5,fp,'--o','linewidth',1.5);
set(gca,'xdir','reverse')
legend({'Empirical precision','Ledoit-Wolf','Permutation Testing','Mixture modeling','Hard threshold, top 5%','Hard threshold, top 10%'},'location','southwest','FontSize',8)
xlabel('Number of nodes in network','FontSize',12)
ylabel('False positive rate','FontSize',12)

f3= figure('visible','off');
plot(49:-1:5,tp,'--o','linewidth',1.5);
set(gca,'xdir','reverse')
legend({'Empirical precision','Ledoit-Wolf','Permutation Testing','Mixture modeling','Hard threshold, top 5%','Hard threshold, top 10%'},'location','southwest','FontSize',8)
xlabel('Number of nodes in network','FontSize',12)
ylabel('True positive rate','FontSize',12)

saveTightFigure(f1,'plots/varN/varNp')
saveTightFigure(f2,'plots/varN/varNf')
saveTightFigure(f3,'plots/varN/varNt')

end