function plotROC
% PLOTROCCURVE Calculates ROC curve for select set of parameters
% Produces: Figure 1D for the paper:
% Thresholding functional connectomes by means of Mixture Modeling (2017)  
% by Bielczyk, Walocha et al.
% Author: Fabian Walocha (2017), fawalocha@gmail.com

clc;
clearvars;
close all;


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
connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
conns = sum(sum(connM));

%% Normalization of data

for IDX = 1:subjs
    tmp = data(IDX,:,:);
    data(IDX,:,:) = (tmp-mean(tmp(:)))./std(tmp(:));
end

%% Calculation inverse covariance matrix

[par_corrs1,par_corrs2] = deal(zeros([subjs,roi,roi]));


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

%% Permutation test

% Permutation testing
maxShuffle = 1000;
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

res_permT2 = zeros([subjs,roi,roi]);
for IDX1 = 1:subjs
    for IDX2 = 1:roi
        for IDX3 = 1:roi
            res_permT2(IDX1,IDX2,IDX3) = ...
                sum(perm_covs(:,IDX2,IDX3)>par_corrs1(IDX1,IDX2,IDX3))/maxShuffle;       
        end
    end
end

paramV = 0:0.1:6;
paramV = paramV/6;
paramV = paramV.^2;
paramV(62) = 0.05;
[tp,fp] = deal(zeros([4,length(paramV)]));

[res_N1,res_N2] = deal(par_corrs1);
for IDX = 1:subjs
    tmp = squeeze(par_corrs1(IDX,:,:));
    tmp = sort(tmp(connM),'descend');
    res_N1(IDX,:,:) = par_corrs1(IDX,:,:) > tmp(floor(conns*0.05));
    res_N2(IDX,:,:) = par_corrs1(IDX,:,:) > tmp(floor(conns*0.1));
end

[tpr5,tpr6,fpr5,fpr6] = deal(zeros(subjs,1));

for IDX = 1:subjs
    tmp5 = squeeze(res_N1(IDX,:,:));
    tmp6 = squeeze(res_N2(IDX,:,:));
    tmpnet = squeeze(net(IDX,:,:));

    tmp5 = tmp5(connM);
    tmp6 = tmp6(connM);
    tmpnet = tmpnet(connM);

    tpr5(IDX) = sum(tmp5.*tmpnet)/sum(tmpnet~=0);
    tpr6(IDX) = sum(tmp6.*tmpnet)/sum(tmpnet~=0);
    
    fpr5(IDX) = sum(tmp5.*(tmpnet==0))/sum(tmpnet==0);
    fpr6(IDX) = sum(tmp6.*(tmpnet==0))/sum(tmpnet==0);
end

tp2(3) = mean(tpr5);
tp2(4) = mean(tpr6);

fp2(3) = mean(fpr5);
fp2(4) = mean(fpr6);

for index = 1:length(paramV)
    %% Other methods
    
    res_MM3 = do_MM2(par_corrs2,paramV(index),'GGM');    

    res_permT = res_permT2<paramV(index);
    res_EP = par_corrs1>sqrt(paramV(index))*2-1;
    res_ledW = par_corrs2>sqrt(paramV(index))*2-1;

    [tpr1,tpr2,tpr3,tpr4,fpr1,fpr2,fpr3,fpr4] = deal(zeros(subjs,1));

    for IDX = 1:subjs
        tmp1 = squeeze(res_permT(IDX,:,:));
        tmp2 = squeeze(res_MM3(IDX,:,:));
        tmp3 = squeeze(res_EP(IDX,:,:));
        tmp4 = squeeze(res_ledW(IDX,:,:));
        tmpnet = squeeze(net(IDX,:,:));

        tmp1 = tmp1(connM);
        tmp2 = tmp2(connM);
        tmp3 = tmp3(connM);
        tmp4 = tmp4(connM);
        tmpnet = tmpnet(connM);

        tpr1(IDX) = sum(tmp1.*tmpnet)/sum(tmpnet~=0);
        tpr2(IDX) = sum(tmp2.*tmpnet)/sum(tmpnet~=0);
        tpr3(IDX) = sum(tmp3.*tmpnet)/sum(tmpnet~=0);
        tpr4(IDX) = sum(tmp4.*tmpnet)/sum(tmpnet~=0);

        fpr1(IDX) = sum(tmp1.*(tmpnet==0))/sum(tmpnet==0);
        fpr2(IDX) = sum(tmp2.*(tmpnet==0))/sum(tmpnet==0);
        fpr3(IDX) = sum(tmp3.*(tmpnet==0))/sum(tmpnet==0);
        fpr4(IDX) = sum(tmp4.*(tmpnet==0))/sum(tmpnet==0);
    end
        
    tp(1,index) = mean(tpr1);
    tp(2,index) = mean(tpr2);
    tp(3,index) = mean(tpr3);
    tp(4,index) = mean(tpr4);
    
    fp(1,index) = mean(fpr1);
    fp(2,index) = mean(fpr2);
    fp(3,index) = mean(fpr3);
    fp(4,index) = mean(fpr4);
end


names = {'Permutation testing','Mixture modeling','Empirical precision','Ledoit-Wolf',...
         'Permutation testing,p<0.05','Mixture modeling,FDR<0.05',...
         'Proportional thresholding, 5%', 'Proportional thresholding, 10%'};

fig1 = figure('visible','off');
hold on;
for IDX = [1 2 3 4]
    a = fp(IDX,1:61);
    b = tp(IDX,1:61);
    plot(a,b,'--.')
end

plot(fp(1,62),tp(1,62),'o','LineWidth',1.5)
plot(fp(2,62),tp(2,62),'o','LineWidth',1.5)
plot(fp2(3),tp2(3),'o','LineWidth',1.5)
plot(fp2(4),tp2(4),'o','LineWidth',1.5)

plot([0,1],[0,1],'--k')
xlim([0,0.2])
xlabel('FPR','Fontsize',12)
ylabel('TPR','Fontsize',12)
legend(names,'location','southeast','fontsize',5)
title('ROC curve for benchmark dataset','Fontsize',10)

saveTightFigure(fig1,'plots/ROCcurve');
end
