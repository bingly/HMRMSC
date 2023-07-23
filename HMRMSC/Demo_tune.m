%% For convinience, we assume the order of the tensor is always 3;
clear;
addpath('tSVD','proxFunctions','solvers','twist');
addpath('ClusteringMeasure', 'LRR', 'Nuclear_norm_l21_Algorithm', 'unlocbox','tools');
addpath(genpath('gspbox'));

dataname='MSRC-v1';  %EYaleB10_mtv scene15 MSRC-v1 BBCSport COIL20_3VIEWS reutersMulSubset UCI_3view ORL_mtv

%% 
%data preparation...
cd ..
[X, gt, K] = data_preparation(dataname);
cd ./HMRMSC

% beta = [1, 10, 100]';  %加权系数
for v = 1:K
    X{v}=NormalizeData(X{v});
    beta(v) = 10^(v-1);
end


% % authentic data in papaer
% lambda1 = 0.4;
% lambda2 = 0.4;
% lambda3 = 0.4;
% 
% JRL_hyper(X, gt, lambda1, lambda2, lambda3, beta);
% 
% JRL_hyper_iter(X, gt, lambda1, lambda2, lambda3, beta, 1, true);


%% iter for best performance
lambda = [0.4,0.8,1.2,1.6,2];
for lambda1 = lambda
    for lambda2 = lambda
        for lambda3 = lambda
            JRL_hyper(X, gt, lambda1, lambda2, lambda3, beta);
        end
    end
end


 

