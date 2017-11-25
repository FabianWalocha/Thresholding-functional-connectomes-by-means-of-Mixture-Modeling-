function plotBIC
% PLOTBIC Comparison of all 8 methods on synthetic, HCP rest, and HCP task via BIC. 
% Produces: Figure 6, figure 7, figure 8 for the paper:
% Thresholding functional connectomes by means of Mixture Modeling (2017)  
% by Bielczyk, Walocha et al.
% Author: Fabian Walocha (2017), fawalocha@gmail.com

clearvars;
close all;
clc;

for IDX10 = 1:3

    if IDX10 == 1
        %% Simulation data

        % Read data
        data_mat = importdata(strcat(pwd,'/data/Smith_orig/sim4.mat'));

        subjs = data_mat.Nsubjects;
        roi = data_mat.Nnodes;
        ts = size(data_mat.ts,1)/subjs;
        net = data_mat.net;
        net = net ~= 0;
        connM = logical(triu(ones(roi))-diag(diag(ones(roi))));
        conns = sum(sum(connM));


        % Format data into SUBJSxTSxROI matrix
        data = zeros([subjs,ts,roi]);
        for IDX1 = 1:subjs
            data(IDX1,:,:) = real(data_mat.ts((1+ts*(IDX1-1)):(ts*IDX1),:));
        end
    elseif IDX10 == 2
        %% Real data (rest)

        cort_sys = 'visual';
        dataset = 'rest';
        
        pwdata = strcat('/data/HCP_',cort_sys,'_highpassfiltered_smoothed3mm_',dataset,'/');
        tmt = csvread('data/SubjsParcUnr.txt');

        data = cell(length(tmt),1);
        for IDX = 1:length(tmt)
            data{IDX} = load(strcat(pwd,pwdata,cort_sys,num2str(tmt(IDX)),'.mat')); 
        end

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

        data_mat_short = data_mat(:,:,:,:,[1:11 14:20]);

        d_short = cat(5,data_mat_short(:,1,:,:,:),data_mat_short(:,2,:,:,:));
        d_short = squeeze(d_short);

        d_short = d_short(:,:,[6:205 1206:1405],:);
        data = squeeze(cat(3,d_short(:,1,:,:),d_short(:,2,:,:)));
    else
       %% Real data (task)
       
        cort_sys = 'visual';
        dataset = 'task';

        pwdata = strcat('/data/HCP_',cort_sys,'_highpassfiltered_smoothed3mm_',dataset,'/');

        tmt = csvread('data/SubjsParcUnr.txt');

        data = cell(length(tmt),1);
        for IDX = 1:length(tmt)
            data{IDX} = load(strcat(pwd,pwdata,cort_sys,num2str(tmt(IDX)),'.mat')); 
        end

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

        data_mat_short = data_mat(:,:,:,:,[1:11 14:20]);

        d_short = cat(5,data_mat_short(:,1,:,:,:),data_mat_short(:,2,:,:,:));
        d_short = squeeze(d_short);
        d_short_tmp = d_short(:,:,6:405,:);
        d_short = d_short(:,:,6:405,:);
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


    %% MM
    bics(1,:) = do_BIC(par_corrs1,'GGM');
    bics(2,:) = do_BIC(par_corrs2,'GGM');
    bics(3,:) = do_BIC(par_corrs1,'GIM');
    bics(4,:) = do_BIC(par_corrs2,'GIM');
    bics(5,:) = do_BIC(par_corrs1,'LGM');
    bics(6,:) = do_BIC(par_corrs2,'LGM');
    bics(7,:) = do_BIC(par_corrs1,'LIM');
    bics(8,:) = do_BIC(par_corrs2,'LIM');



    names = {'MM(EP,GG)','MM(LW,GG)','MM(EP,GIG)','MM(LW,GIG)',...
          'MM(EP,LG)','MM(LW,LG)','MM(EP,LIG)','MM(LW,LIG)'};

    fig1=figure('visible','off');
    hold on;
    c = zeros([size(bics,1),3]);
    for IDX = 1:size(bics,1)
        p = plot(bics(IDX,:)','o','linewidth',1.5);
        c(IDX,:) = get(p,'color');
    end
    for IDX = 1:size(bics,1)
        plot([0,subjs],[nanmean(bics(IDX,:)),nanmean(bics(IDX,:))],'color',c(IDX,:),'linewidth',2.5);
    end
    xlabel('Subjects','Fontsize',12)
    legend(names,'Location','southeast')
    ylabel('BIC value')
    xlim([1,subjs])

    a = nanmean(bics([1 3 5 7],:),1);
    b = nanmean(bics([2 4 6 8],:),1);

    vals(1,1) = nanmean(a,2);
    vals(2,1) = nanmean(b,2);
    vals(3,1) = ranksum(a,b);

    fig2= figure('visible','off');
    hold on;
    c = zeros([2,3]);
    p = plot(a,'o','linewidth',1.5);  
    c(1,:) = get(p,'color');
    p = plot(b,'o','linewidth',1.5);
    c(2,:) = get(p,'color');
    plot([0,subjs],[nanmean(a),nanmean(a)],'color',c(1,:),'linewidth',2.5);
    plot([0,subjs],[nanmean(b),nanmean(b)],'color',c(2,:),'linewidth',2.5);
    xlabel('Subjects','Fontsize',12)
    ylabel('BIC value')
    xlim([1,subjs])
    legend('Empirical precision based','Ledoit-Wolf based','Location','southwest')

    a = nanmean(bics(1:4,:),1);
    b = nanmean(bics(5:8,:),1);

    vals(1,2) = nanmean(a,2);
    vals(2,2) = nanmean(b,2);
    vals(3,2) = ranksum(a,b);

    fig3= figure('visible','off');
    hold on;
    c = zeros([2,3]);
    p = plot(a,'o','linewidth',1.5);  
    c(1,:) = get(p,'color');
    p = plot(b,'o','linewidth',1.5);
    c(2,:) = get(p,'color');
    plot([0,subjs],[nanmean(a),nanmean(a)],'color',c(1,:),'linewidth',2.5);
    plot([0,subjs],[nanmean(b),nanmean(b)],'color',c(2,:),'linewidth',2.5);
    xlabel('Subjects','Fontsize',12)
    ylabel('BIC values')
    xlim([1,subjs])
    legend('Gaussian pseudo-null','Laplacian pseudo-null','Location','southwest')

    a = nanmean(bics([1 2 5 6],:),1);
    b = nanmean(bics([3 4 7 8],:),1);

    vals(1,3) = nanmean(a,2);
    vals(2,3) = nanmean(b,2);
    vals(3,3) = ranksum(a,b);

    fig4= figure('visible','off');
    hold on;
    c = zeros([2,3]);
    p = plot(a,'o','linewidth',1.5);  
    c(1,:) = get(p,'color');
    p = plot(b,'o','linewidth',1.5);
    c(2,:) = get(p,'color');
    plot([0,subjs],[nanmean(a),nanmean(a)],'color',c(1,:),'linewidth',2.5);
    plot([0,subjs],[nanmean(b),nanmean(b)],'color',c(2,:),'linewidth',2.5);
    xlabel('Subjects','Fontsize',12)
    ylabel('BIC values')
    xlim([1,subjs])
    legend('Gamma signal','Inverse-Gamma signal','Location','southwest')

    vals(1,4) = nanmean(bics(1,:),2);
    vals(2,4) = nanmean(bics(2,:),2);
    vals(3,4) = nanmean(bics(3,:),2);
    vals(4,4) = nanmean(bics(4,:),2);
    vals(5,4) = nanmean(bics(5,:),2);
    vals(6,4) = nanmean(bics(6,:),2);
    vals(7,4) = nanmean(bics(7,:),2);
    vals(8,4) = nanmean(bics(8,:),2);


    if IDX10 == 1

        saveTightFigure(fig1,'plots/bic/BIC11');
        saveTightFigure(fig2,'plots/bic/BIC21');
        saveTightFigure(fig3,'plots/bic/BIC31');
        saveTightFigure(fig4,'plots/bic/BIC41');
        save('results/bicVals1','vals');
    elseif IDX10 == 2
        saveTightFigure(fig1,'plots/bic/BIC12');
        saveTightFigure(fig2,'plots/bic/BIC22');
        saveTightFigure(fig3,'plots/bic/BIC32');
        saveTightFigure(fig4,'plots/bic/BIC42');
        save('results/bicVals2','vals');
    else
        saveTightFigure(fig1,'plots/bic/BIC13');
        saveTightFigure(fig2,'plots/bic/BIC23');
        saveTightFigure(fig3,'plots/bic/BIC33');
        saveTightFigure(fig4,'plots/bic/BIC43');
        save('results/bicVals3','vals');
    end
    clearvars bics
end
end

