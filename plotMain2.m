function [dens] = plotMain2( dataRM, dataTM,dataRP,dataTP,dataRE,dataTE,dataRL,dataTL,dataRN1,dataTN1,dataRN2,dataTN2)
%PLOTMAIN2 second part of plotting, used for getting difference matrices
%   Plots heatmaps for permutation test and mixture modelling in task and
%   rest (and difference)
%   Inputs:
%       dataRM  rest mixture modeling
%       dataTM  task mixture modeling
%       dataRP  rest permutation testing
%       dataTP  task permutation testing
%       dataRE  rest empirical precision
%       dataTE  task empirical precision
%       dataRL  rest ledoit-wolf
%       dataTL  task ledoit-wolf
%       dataRN1 rest proportional thresholding 5%
%       dataTN1 task proportional thresholding 5%
%       dataRN2 rest proportional thresholding 10%
%       dataTN2 task proportional thresholding 10%
%   Outputs:
%       dens    density of differences and significant differences
%
% Saves codes for figure 4, figure 5, figure 11
% Author: Fabian Walocha (2017), fawalocha@gmail.com

nams = {'EP','ledW','permT','MM','N1','N2'};

dens = zeros([6,2]);
for IDX = 1:6
    switch IDX
        case 1
            resultR = dataRE;
            resultT = dataTE;
        case 2
            resultR = dataRL;
            resultT = dataTL;
        case 3
            resultR = dataRP;
            resultT = dataTP;
        case 4 
            resultR = dataRM;
            resultT = dataTM;
        case 5
            resultR = dataRN1;
            resultT = dataTN1;
        case 6
            resultR = dataRN2;
            resultT = dataTN2;
    end
    
    resultD = squeeze(sum(resultT-resultR,1));
    mwu2 = mannwhitu(resultR,resultT);
    resultD3 = resultD.*mwu2;
    
    dens(IDX,1) = sum(sum(triu(resultD~=0)))/630;
    dens(IDX,2) = sum(sum(triu(mwu2~=0)))/630;
    
    switch IDX
        case 1
            res_EP = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_EP');

            res_EP = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_EP');
            
        case 2
            res_ledW = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_ledW');

            res_ledW = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_ledW');
        case 3
            res_permT = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_permT');

            res_permT = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_permT');
        case 4
            res_MM = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_MM');

            res_MM = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_MM');
        case 5
            res_N1 = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_N1');

            res_N1 = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_N1');
        case 6
            res_N2 = resultD;
            wx = strcat('results/res_',nams{IDX},'diff');
            save(wx,'res_N2');

            res_N2 = resultD3;
            wx = strcat('results/res_',nams{IDX},'diff3');
            save(wx,'res_N2');
    end

end
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