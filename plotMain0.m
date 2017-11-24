function plotMain0
%PLOTMAIN0 Arranges data in order to prepare for plotMain1
% Author: Fabian Walocha (2017), fawalocha@gmail.com

%% Prepare data

clc;
clearvars;
close all;

%% Set parameters

% real data or simulation data
dat = 'real';
% dat = 'sim';
mmtype = 'GGM';
lwFlag = 1;

% cutoff fdr value 
fdrV = 0.05;

if strcmp(dat,'sim')

    %% Simulation data

    % Read data
    data_mat = importdata(strcat(pwd,'/data/Smith_orig/sim4.mat'));

    subjs = data_mat.Nsubjects;
    roi = data_mat.Nnodes;
    ts = size(data_mat.ts,1)/subjs;
    net = data_mat.net;
    net = net ~= 0;


    % Format data into SUBJSxTSxROI matrix
    data = zeros([subjs,ts,roi]);
    for IDX1 = 1:subjs
        data(IDX1,:,:) = real(data_mat.ts((1+ts*(IDX1-1)):(ts*IDX1),:));
    end
    plotMain1(data,fdrV,lwFlag,0,mmtype,net);

elseif strcmp(dat,'real')
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

    [rRM,rRP,rRE,rRL,rRN1,rRN2,iR,avgsR,i2R] = plotMain1(data,fdrV,lwFlag,1,mmtype);
    
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

    [rTM,rTP,rTE,rTL,rTN1,rTN2,iT,avgsT,i2T] = plotMain1(data,fdrV,lwFlag,2,mmtype);
    
    %% Save further data on results
    
    % Average number of connections per subject
    avgs = [avgsR',avgsT'];
    
    standev = zeros(6,2);
    for IDX = 1:6
        standev(IDX,1) = std(iR(IDX,:));
        standev(IDX,2) = std(iT(IDX,:));
    end
    
    % index overlap and ICC  scores
    iccs = [mean(iR,2),mean(iT,2),standev];
    iccs2 = [i2R,i2T];
    
    %% 

    % dens = density of group connectomes 
    dens = plotMain2(rRM,rTM,rRP,rTP,rRE,rTE,rRL,rTL,rRN1,rTN1,rRN2,rTN2);
    save('results/iccs','iccs');
    save('results/avgs','avgs');
    save('results/dens','dens');
    save('results/iccs2','iccs2');
    
end
end