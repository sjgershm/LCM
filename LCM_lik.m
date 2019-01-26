function [lik, latents] = LCM_lik(alpha,data,opts)
    
    % Compute log-likelihood of data under the latent cause model.
    %
    % USAGE: [lik, latents] = LCM_lik(alpha,data,[opts])
    %
    % INPUTS:
    %   alpha - concentration parameter
    %   data - single-subject data
    %   opts - options structure
    %
    % OUTPUTS:
    %   lik - log-likelihood
    %   latents - structure containing latent variables:
    %               .b - beta coefficient mapping model CR to  measured CR
    %               .sd - maximum likelihood standard deviation
    %               .CR - predicted CR
    %
    % Sam Gershman, Jan 2019
    
    % set concentration parameter
    opts.alpha = alpha;
    
    % run particle filter
    results = LCM_infer([data.US data.CS],opts);
    
    % use linear regression to fit model output to CR
    N = length(results.V);
    X = results.V;
    b = (X'*X)\(X'*data.CR);                % maximum likelihood regression coefficients
    CRpred = X*b;                           % predicted CR
    sd = sqrt(mean((data.CR - CRpred).^2)); % maximum likelihood standard deviation
    lik = sum(log(normpdf(data.CR,CRpred,sd)));  % log-likelihood
    
    % return latent variables
    if nargout > 1
        latents = results;
        latents.b = b;
        latents.sd = sd;
        latents.CR = CRpred;
    end