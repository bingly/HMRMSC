function JRL_hyper_iter(X, gt, lambda1,lambda2,lambda3, beta, count, cal_time)
% JRL_hyper_iter
% Authors: Bingcai 10/15/2021
% beta 权重奇异值
% count 总共运行次数--计算平均值和标准差
% cal_time 调试打印时间

%% Note: each column is an sample (same as in LRR)
%% 

% Initialize...
K = length(X); N = size(X{1},2); %sample number

%% 循环开始
for index=1:count
    if(exist('cal_time', 'var'))
        tic;
    end
    
    for k=1:K
        Z{k} = zeros(N,N); %Z{2} = zeros(N,N);
        L{k} = zeros(N,N);
        S{k} = zeros(N,N); %G
        E{k} = zeros(size(X{k},1),N); %E{2} = zeros(size(X{k},1),N);
        Y{k} = zeros(size(X{k},1),N); %Y{2} = zeros(size(X{k},1),N);
        XXv{k} = X{k}' * X{k};
    end

%     z = zeros(N*N*K,1);
%     s = zeros(N*N*K,1);

    sX = [N, N, K];
    %set Default
    parOP         =    false;
    ABSTOL        =    1e-6;
    RELTOL        =    1e-4;


    Isconverg = 0;epson = 1e-7;
    iter = 0;
    mu = 1; max_mu = 10e10; pho_mu = 2;
    start = 1;

    while(Isconverg == 0)
    %     fprintf('----processing iter %d--------\n', iter+1);
    %-------------------0 update L^k-------------------------------
        for k=1:K
            if start==1
              Weight{k} = constructW_PKN((abs(Z{k})+abs(Z{k}'))./2, 3);
              Diag_tmp = diag(sum(Weight{k}));
              L{k} = Diag_tmp - Weight{k};
            else
            %------------modified to hyper-graph---------------
              P =  (abs(Z{k})+abs(Z{k}'))./2;
              param.k = 3;
              HG = gsp_nn_hypergraph(P', param);
              L{k} = HG.L;
            end

    %         Weight{k} = constructW_PKN((abs(Z{k})+abs(Z{k}'))./2, 10);
    %         Diag_tmp = diag(sum(Weight{k}));
    %         L{k} = Diag_tmp - Weight{k};
        end
        start = 0;

        for k=1:K
            %1 update Z^k
            A = lambda3*eye(N) + mu*XXv{k};
            B = 2*lambda1*L{k};
            C = lambda3*S{k} + X{k}'*Y{k} + mu*XXv{k} - mu*X{k}'*E{k};
            Z{k} = lyap(A,B,-C);

            %2 update E^k
            %origin
            %F = [X{1}-X{1}*Z{1}+Y{1}/mu;X{2}-X{2}*Z{2}+Y{2}/mu;X{3}-X{3}*Z{3}+Y{3}/mu];
            
            %use iteral
            %first cal K>=2, then cal K>=3
            F = [X{1}-X{1}*Z{1}+Y{1}/mu;X{2}-X{2}*Z{2}+Y{2}/mu];
            for k1=3:K
                F = [F;X{k1}-X{k1}*Z{k1}+Y{k1}/mu];
            end
            
            %[Econcat,info] = prox_l21(F, 0.5/1);
            [Econcat] = solve_l1l2(F,1/mu);
            e_start = 1; e_end = 0;
            for k1 = 1:K
                e_end = e_end + size(X{k1},1);
                E{k1} = Econcat(e_start:e_end, :);
                e_start = e_start + size(X{k1},1);
            end
%             E{1} = Econcat(1:size(X{1},1),:);
%             E{2} = Econcat(size(X{1},1)+1:size(X{1},1)+size(X{2},1),:);
%             E{3} = Econcat(size(X{1},1)+size(X{2},1)+1:end,:);

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
                err(k) = history.norm_Z;
    %             fprintf('RR    norm_Z %7.10f    ', history.norm_Z);
                Isconverg = 0;
            end

            S{k} = S_tensor(:,:,k);
        end
        fprintf('iter-%2d: %7.10f\n',iter+1 ,mean(err));

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
    
%     figure(1); imagesc(Aff_S); %affinity matrix

     cls_num = length(unique(gt));
     C = SpectralClustering(Aff_S,cls_num);
     
    [A nmi(index) avgent] = compute_nmi(gt,C);
    ACC(index) = Accuracy(C,double(gt));
    [f(index),p(index),r(index)] = compute_f(gt,C);
    [AR(index),RI,MI,HI]=RandIndex(gt,C);
    
    % fprintf('\n%.4f %.4f %.4f %.4f %.4f %.4f\n',ACC,nmi,AR,f,p,r);
    fprintf('%d-%.2f,%0.2f,%.2f: %.4f %.4f %.4f %.4f %.4f %.4f\n',index,lambda1,lambda2,lambda3,ACC(index),nmi(index),AR(index),f(index),p(index),r(index));
    if(exist('cal_time', 'var'))
        toc;
    end
end

fprintf('HMRMSC,%.3f±%.3f,%.3f±%.3f,%.3f±%.3f,%.3f±%.3f,%.3f±%.3f,%.3f±%.3f \n',...
    mean(ACC),std(ACC),mean(nmi),std(nmi),mean(AR),std(AR),mean(f),std(f),mean(p),std(p),mean(r),std(r));
