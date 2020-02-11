function results = LCM_infer(X,opts)
    
    % Particle filtering or local maximum a posteriori inference for latent cause
    % model of associative learning.
    %
    % USAGE: results = LCM(X,[opts])
    %
    % INPUTS:
    %   X - [T x D] stimulus inputs, where T is the number of timepoints
    %       and D is the number of stimulus features. We assume binary
    %       features. The first feature (column 1) is the US, and the rest
    %       of the features (column 2 through D) are CSs.
    %   opts (optional) - structure containing various options (see
    %       LCM_opts). Missing fields are set to defaults. If opts.M = 1,
    %       then the model computes a local MAP estimate.
    %
    % OUTPUTS:
    %   results - structure wih the following fields:
    %       .opts - options (missing fields set to defaults)
    %       .V - [T x 1] US prediction
    %       .post - [T x K] latent cause posterior, where post(t,k) is the
    %            probability of latent cause k being active on trial t,
    %            after observing the all the features. K (the number of
    %            latent causes) is determined adaptively by the model.
    %
    % Sam Gershman, July 2016
    
    % set parameters
    if nargin < 2; opts = []; end
    opts = LCM_opts(opts);
    M = opts.M;
    a = opts.a;
    b = opts.b;
    results.opts = opts;
    
    % initialization
    if opts.alpha==0; K = 1; else K = opts.K; end
    post = zeros(1,K); post(1) = 1;
    post0 = zeros(M,K); post0(:,1) = 1;
    [T, D] = size(X);
    N = zeros(M,K,D);                           % feature-cause co-occurence counts
    B = zeros(M,K,D);                           % feature-cause co-occurence counts
    Nk = zeros(M,K);                            % cause counts
    results.post = [ones(T,1) zeros(T,K-1)];    % cause assignments
    results.V = zeros(T,1);                     % US predictions
    z = ones(M,1);
    
    % loop over trials
    for t = 1:T
        
        % calculate likelihood
        lik = N;
        lik(:,:,X(t,:)==0) = B(:,:,X(t,:)==0);
        lik = bsxfun(@rdivide,lik+a,Nk+a+b);
        
        if opts.alpha > 0    % only update posterior if concentration parameter is non-zero
            
            % calculate CRP prior
            prior = Nk;
            for m = 1:M
                prior(m,z(m)) = prior(m,z(m)) + opts.stickiness; % add stickiness
                prior(m,find(prior(m,:)==0,1)) = opts.alpha;     % probability of a new latent cause
            end
            prior = bsxfun(@rdivide,prior,sum(prior,2));
            
            % posterior conditional on CS only
            post = prior.*squeeze(prod(lik(:,:,2:D),3));
            post0 = bsxfun(@rdivide,post,sum(post,2));
            
            % posterior conditional on CS and US
            post = post.*squeeze(lik(:,:,1));
            post = post / sum(sum(post));
        end
        results.post(t,:) = mean(bsxfun(@rdivide,post,sum(post,2)),1);
        
        % posterior predictive mean for US
        pUS = squeeze(N(:,:,1)+a)./(Nk+a+b);
        results.V(t,1) = post0(:)'*pUS(:)./M;
        
        % sample new particles
        x1 = X(t,:)==1; x0 = X(t,:)==0;

        if M==1
            [~,z] = max(post);                         % maximum a posteriori
            Nk(1,z) = Nk(1,z) + 1;
            N(1,z,x1) = N(1,z,x1) + 1;
            B(1,z,x0) = B(1,z,x0) + 1;
        else
            NkOld = Nk;
            NOld = N;
            BOld = B;
            for m = 1:M
                row = min(find(rand() < cumsum(sum(post,2))));
                Nk(m,:) = NkOld(row,:);
                N(m,:,:) = NOld(row,:,:);
                B(m,:,:) = BOld(row,:,:);
                col = min(find(rand() < cumsum(post(row,:)/sum(post(row,:)))));
                Nk(m,col) = Nk(m,col) + 1;
                N(m,col,x1) = N(m,col,x1) + 1;
                B(m,col,x0) = B(m,col,x0) + 1;
            end
        end
    end
    
    % remove unused particles
    ix = mean(results.post)==0;
    results.post(:,ix) = [];