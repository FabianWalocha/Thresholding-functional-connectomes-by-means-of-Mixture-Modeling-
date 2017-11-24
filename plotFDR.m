function plotFDR
% Plots the estimated FDR vs the true underlying FDR on synthetic data
% Produces: Figure 1C
% Author: Fabian Walocha (2017), fawalocha@gmail.com

%% Prepare data
clc;
clearvars;
close all;
%% Set parameters
mmtype = 'GGM';
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
par_corrs2 = zeros([subjs,roi,roi]);

    
% Precision values to partial correlations
for IDX1 = 1:subjs
    [sigma_hat,~] = covCor(squeeze(data(IDX1,:,:))); 
    par_corrs2(IDX1,:,:) = inv(sigma_hat);
    for IDX2 = 1:roi
        for IDX3 = IDX2+1:roi
            par_corrs2(IDX1,IDX2,IDX3) = -par_corrs2(IDX1,IDX2,IDX3) / ...
                sqrt(par_corrs2(IDX1,IDX2,IDX2)*par_corrs2(IDX1,IDX3,IDX3));
            par_corrs2(IDX1,IDX3,IDX2) = -par_corrs2(IDX1,IDX3,IDX2) / ...
                sqrt(par_corrs2(IDX1,IDX2,IDX2)*par_corrs2(IDX1,IDX3,IDX3));
        end
    end
end


numP = 100;
fdrT = zeros([subjs,numP]);
paramV = linspace(0,1,numP);
for IDX10 = 1:numP
    if mod(IDX10,10)==0 
        disp(IDX10)
    end
    res_MM = do_MM2(par_corrs2,paramV(IDX10),mmtype);
    for IDX = 1:subjs
        pplus = sum(res_MM(IDX,connM))/conns;
        pN0 = sum(net(IDX,connM)==0)/conns;
        fprT = sum((net(IDX,connM)==0).*(res_MM(IDX,connM)~=0)) / sum(net(IDX,connM)==0);
        fdrT(IDX,IDX10) = (fprT*pN0)/pplus;
    end
end

fdr1 = mean(fdrT,1);
save('results/fdrDats','fdr1');
f1 = figure('visible','off');
hold on;
plot(paramV,paramV,'--','linewidth',1.5)
plot(paramV,fdr1,'linewidth',2)
xlabel('True FDR');
ylabel('Estimated pFDR');
saveTightFigure(f1,'plots/fdr');

end
