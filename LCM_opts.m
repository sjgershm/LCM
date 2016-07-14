function opts = LCM_opts(opts)
    
    % Set options.
    %
    % USAGE: opts = LCM_opts([opts])
    %
    % INPUTS:
    %   opts (optional) - options structure with a subset of fields
    %       specified. All missing or empty fields will be set to defaults. If
    %       opts = [], then the entire structure is set to defaults.
    %
    % OUTPUTS:
    %   opts - fully specified options structure
    %
    % DEFAULTS:
    %   opts.M = 100        (number of particles)
    %   opts.a = 1          (hyperparameter of beta prior: pseudo-count for feature presence)
    %   opts.b = 1          (hyperparameter of beta prior: pseudo-count for feature absence)
    %   opts.alpha = 0      (concentration parameter for Chinese restaurant process prior)
    %   opts.stickiness = 0 (stickiness parameer for Chinese restaurant process prior)
    %   opts.K = 10         (maximum number of latent causes)
    %
    % Sam Gershman, July 2016
    
    def_opts.M = 100;
    def_opts.a = 1;
    def_opts.b = 1;
    def_opts.alpha = 0;
    def_opts.stickiness = 0;
    def_opts.K = 10;
    
    if nargin < 1 || isempty(opts)
        opts = def_opts;
    else
        F = fieldnames(def_opts);
        for f = 1:length(F)
            if ~isfield(opts,F{f}) || isempty(opts.(F{f}))
                opts.(F{f}) = def_opts.(F{f});
            end
        end
    end
    
    % make sure parameters aren't negative
    opts.a = max(opts.a,0);
    opts.b = max(opts.b,0);
    opts.alpha = max(opts.alpha,0);
    opts.stickiness = max(opts.stickiness,0);