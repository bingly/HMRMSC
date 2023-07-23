function [X, gt, K] = data_preparation(dataname)

load(strcat('./data/',dataname,'.mat'));
fprintf('start: %s\n',dataname);

%% Note: each column is an sample (same as in LRR)
%% 
%data preparation...
% X
if(strcmp(dataname, 'UCI_3view'))
    X = data;
    gt = truth;
    clear data;
    clear truth;
end
if(strcmp(dataname, 'COIL20_3VIEWS'))
   X{1}=X1;X{2}=X2;X{3}=X3;X{4}=X4;
end
if(strcmp(dataname, 'scene15'))
   X{1}=X1;X{2}=X2;X{3}=X3;
end
% BBCSport data-preparation
if(exist('fea', 'var'))
    X = fea;
    clear fea;
end

K = length(X);

 % compatible for tag X, ajust to columns represent samples
 if(size(X{1},2) ~= size(X{2},2))
    for k=1:K
        X{k} = X{k}';
    end
end
 
% compatible for tag Y
if(~exist('gt','var'))
    gt = Y;
    clear Y;
end
if(ismember(0,gt))
    gt = gt + 1;
end