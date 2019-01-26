function results = LCM_fit(data,opts)
    
    % Fit latent cause model to data.
    %
    % USAGE: results = LCM_fit(data,[opts])
    %
    % INPUTS:
    %   data - [nSubjects x 1] structure containing the following fields:
    %           .CR - [nTrials x 1] conditioned response
    %           .CS - [nTrials x nCues] conditioned stimului
    %           .US - [nTrials x 1] unconditioned response
    %   opts (optional) - structure defining LCM options (see LCM_opts.m)
    %
    % OUTPUTS:
    %   results - [nSubjects x 1] structure containing the following fields:
    %               .alpha - concentration parameter values
    %               .P - posterior probability distribution over alpha
    %               .lik - log-likelihood for each alpha value
    %               .latents - latent variables for each alpha value (see LCM_infer)
    %               .logBF - log Bayes factor for the alpha>=0 model
    %                       relative to the alpha=0 model
    %
    % Sam Gershman, Jan 2019
    
    N = 50; % number of alpha values to evaluate
    alpha = linspace(0,10,N);
    if nargin < 2; opts = []; end
    
    for s = 1:length(data)
        disp(['Subject ',num2str(s)]);
        for i = 1:length(alpha)
            [results(s).lik(i,1), results(s).latents(i)] = LCM_lik(alpha(i),data(s),opts);
        end
        L = logsumexp(results(s).lik);
        results(s).alpha = alpha;
        results(s).P = exp(results(s).lik-L);
        results(s).alpha = alpha*results(s).P;
        results(s).logBF = L - log(N) - results(s).lik(1);
    end