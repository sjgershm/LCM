function [results, bms_results] = LCM_fit(data)
    
    % Fit variants of the latent cause model to data.
    % Model 1: multiple latent causes with stickiness
    % Model 2: multiple latent causes without stickiness
    % Model 3: single latent cause
    %
    % USAGE: [results, bms_results] = LCM_fit(data)
    %
    % INPUTS:
    %   data - [nSubjects x 1] structure containing the following fields:
    %           .CR - [nTrials x 1] conditioned response
    %           .CS - [nTrials x nCues] conditioned stimului
    %           .US - [nTrials x 1] unconditioned response
    %
    % OUTPUTS:
    %   results - [nModels x 1] structure containing results of fitted models (see
    %           mfit_optimize.m for more info). Also contains a 'latents' subfield
    %           containing latent variables.
    %   bms_results - Bayesian model selection results (see mfit_bms.m for
    %           more info)
    %
    % Sam Gershman, July 2016
    
    % parameter info
    param(1).name = 'alpha';    % concentration parameter
    param(1).logpdf = @(x) 0;   % uniorm prior
    param(1).lb = 0;            % lower bound
    param(1).ub = 10;           % upper bound
    
    param(2).name = 'a';        % beta distribution feature presence pseudo-count
    param(2).logpdf = @(x) 0;   % uniorm prior
    param(2).lb = 0.5;          % lower bound
    param(2).ub = 1.5;          % upper bound
    
    param(3).name = 'stickiness'; % stickiness of chinese restaurant process
    param(3).logpdf = @(x) 0;   % uniorm prior
    param(3).lb = 0;            % lower bound
    param(3).ub = 5;            % upper bound
    
    % model parameterizations (which parameters to include in each model)
    models{1} = 1:3;
    models{2} = 1:2;
    models{3} = 2;
    
    % model options
    opts = LCM_opts;
    
    for j = 1:length(models)
        
        % fit model
        f = @(x,data) LCM_lik(x,data,param(models{j}),opts);    % log-likelihood function
        r = mfit_optimize(f,param(models{j}),data);
        
        % get latent variables
        for i = 1:length(data)
            [~,r.latents(i)] = f(r.x(i,:),data(i));
        end
        
        results(j) = r;
    end
    
    % Bayesian model selection
    bms_results = mfit_bms(results);