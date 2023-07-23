function JRL(X, gt, lambda1,lambda2,lambda3, beta,cal_time)

if(exist('cal_time', 'var'))
    tic;
end
%% Note: each column is an sample (same as in LRR)
%% 
%data preparation...
 for v=1:3
    [X{v}]=NormalizeData(X{v});
     %X{v} = zscore(X{v},1);
end
% Initialize...

K = length(X); N = size(X{1},2); %sample number

for k=1:K
    Z{k} = zeros(N,N); %Z{2} = zeros(N,N);
    L{k} = zeros(N,N);
    S{k} = zeros(N,N); %G
    E{k} = zeros(size(X{k},1),N); %E{2} = zeros(size(X{k},1),N);
    Y{k} = zeros(size(X{k},1),N); %Y{2} = zeros(size(X{k},1),N);
    XXv{k} = X{k}' * X{k};
end

w = zeros(N*N*K,1);
g = zeros(N*N*K,1);
dim1 = N;dim2 = N;dim3 = K;
myNorm = 'tSVD_1';
sX = [N, N, K];
%set Default
parOP         =    false;
ABSTOL        =    1e-6;
RELTOL        =    1e-4;

%construct L
options = [];
options.WeightMode = 'HeatKernel';  
options.NeighborMode = 'KNN';
options.k = 5;
options.t = 10;
for k=1:K
    W = constructW(X{k}', options);
    DCol = full(sum(W,2));
    D = spdiags(DCol,0,N,N);
    L{k} = D - W;
end

Isconverg = 0;epson = 1e-7;
iter = 0;
mu = 1; max_mu = 10e10; pho_mu = 2;


while(Isconverg == 0)
%     fprintf('----processing iter %d--------\n', iter+1);
    for k=1:K
        %1 update Z^k
        A = lambda3*eye(N) + mu*XXv{k};
        B = 2*lambda1*L{k};
        C = lambda3*S{k} + X{k}'*Y{k} + mu*XXv{k} - mu*X{k}'*E{k};
        Z{k} = lyap(A,B,-C);
        
        %2 update E^k
        F = [X{1}-X{1}*Z{1}+Y{1}/mu;X{2}-X{2}*Z{2}+Y{2}/mu;X{3}-X{3}*Z{3}+Y{3}/mu];

        [Econcat] = solve_l1l2(F,1/mu);
        %F = F';
        %[Econcat,info] = prox_l21(F, 0.5/1);
        E{1} = Econcat(1:size(X{1},1),:);
        E{2} = Econcat(size(X{1},1)+1:size(X{1},1)+size(X{2},1),:);
        E{3} = Econcat(size(X{1},1)+size(X{2},1)+1:end,:);
        %3 update Yk
        %Y{k} = Y{k} + mu*(X{k}-X{k}*Z{k}-E{k});
        Y{k} = Y{k} + mu*(X{k}-X{k}*Z{k}-E{k});
    end
    
    %4 update S
    Z_tensor = cat(3, Z{:,:});
    z = Z_tensor(:);
    
    %twist-version
%    [s, objV] = wshrinkObj(z,lambda2/lambda3,sX,0,3);
   [s, objV] = wshrinkObj_weight(z,beta*lambda2/lambda3,sX,0,3);
    S_tensor = reshape(s, sX);
    
    %record the iteration information
    history.objval(iter+1)   =  objV;

    %% coverge condition
    Isconverg = 1;
    for k=1:K
        if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
            history.norm_Z = norm(X{k}-X{k}*Z{k}-E{k},inf);
%             fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end
        
        S{k} = S_tensor(:,:,k);
    end
   
    if (iter>50)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
end

Aff_S = 0;
for k=1:K
    Aff_S = Aff_S + abs(S{k})+abs(S{k}');
end
Aff_S = Aff_S - diag(diag(Aff_S));
Aff_S = Aff_S / (2*length(X));

 cls_num = length(unique(gt));
 C = SpectralClustering(Aff_S,cls_num);
[A nmi avgent] = compute_nmi(gt,C);
ACC = Accuracy(C,double(gt));
[f,p,r] = compute_f(gt,C);
[AR,RI,MI,HI]=RandIndex(gt,C);
if(exist('cal_time', 'var'))
    toc;
end
% fprintf('\n%.4f %.4f %.4f %.4f %.4f %.4f\n',ACC,nmi,AR,f,p,r);
fprintf('\n%.2f,%.2f,%.2f: %.4f %.4f %.4f %.4f %.4f %.4f\n',lambda1,lambda2,lambda3,ACC,nmi,AR,f,p,r);