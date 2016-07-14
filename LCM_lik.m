function [lik, latents] = LCM_lik(x,data,param,opts)
    
    % Compute log-likelihood of data under the latent cause model.
    %
    % USAGE: [lik, latents] = LCM_lik(x,data,param,[opts])
    %
    % INPUTS:
    %   x - parameter vector
    %   data - single-subject data
    %   param - parameter structure
    %   opts - options structure
    %
    % OUTPUTS:
    %   lik - log-likelihood
    %   latents - structure containing latent variables
    %
    % Sam Gershman, July 2016
    
    % set parameters
    for i = 1:length(param)
        opts.(param(i).name) = x(i);
    end
    
    % run particle filter
    results = LCM_infer([data.US data.CS],opts);
    
    % use linear regression to fit model output to CR
    N = length(results.V);
    X = [ones(N,1) results.V];
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