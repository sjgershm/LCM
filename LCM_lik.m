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
    %               .ix - data points used in the analysis (excluding NaNs)
    %
    % Sam Gershman, Jan 2019
    
    % set concentration parameter
    opts.alpha = alpha;
    
    % run particle filter
    results = LCM_infer([data.US data.CS],opts);
    
    % use linear regression to fit model output to CR
    ix = ~isnan(data.CR);	% exclude missing data
    X = results.V(ix);
    b = (X'*X)\(X'*data.CR(ix));                % maximum likelihood regression coefficients
    CRpred = results.V*b;                           % predicted CR
    sd = sqrt(mean((data.CR(ix) - CRpred(ix)).^2)); % maximum likelihood standard deviation
    lik = sum(log(normpdf(data.CR(ix),CRpred(ix),sd)));  % log-likelihood
    
    % return latent variables
    if nargout > 1
        latents = results;
        latents.b = b;
        latents.sd = sd;
        latents.CR = CRpred;
        latents.ix = ix;
    end