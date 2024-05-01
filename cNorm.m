%===============================================================================
%===============================================================================
%
% This function normalizes a time-dependent concentration vector.
%
%===============================================================================
%
% Inputs:
%   t        := time vector
%   c        := vector of concentrations (mg/L)
%   mymeas   := measurement (conservative or reactive)
%   varargin := optional inputs, based on available data and the
%               specific measurement (conservative or reactive)
%               (must be specified as " 'name', value " pairs)
%   Optional:
%       Q     := Discharge. If zero, discharge will be calculated from dilution
%                gauging, assuming full mass recovery. Only add if discharge
%                has been measured independently.
%       cMass := Input mass of conservative tracer (mg)
%       rMass := Input mass of reactive tracer (mg).
%       fmc   := Fraction mass recovery of conservative tracer.
%
% Outputs:
%   cnorm := normalized concentration
%   fm    := fraction of mass recovered (optional)
%   Q     := discharge, estimated from dilution gauging (optional)
%
%   NOTE: cnorm is a required output, while fm and Q will only be output
%           if the function is called with an output vector of length 2 or 3,
%           respectively.
%
%===============================================================================
%===============================================================================

function [cnorm, varargout] = cNorm(t, c, mymeas, varargin)

%% Function Setup
p = inputParser;

% Error checking functions
check_t = @(x) validateattributes(x, {'double'},...
                                  {'increasing', 'nonempty', 'nonnegative'});
check_c = @(x)  validateattributes(x, {'double'}, {'nonempty', 'nonnegative'});
checkMeas = @(x) any(validatestring(x, {'c', 'r'}));


% Required input parameters
addRequired(p, 't', check_t); % time vector
addRequired(p, 'c', check_c); % concentration vector (mg/L)
addRequired(p, 'mymeas', checkMeas); % measurement, must be 'c' (conservative) or 'r' (reactive)

% Optional inputs.
% Discharge. If zero, discharge will be calculated from dilution gauging,
% assuming full mass recovery. Only add if discharge has been measured
% independently.
addParameter(p, 'Q', -inf, @(x) isnumeric(x) && isscalar(x));

% Input mass of conservative tracer (mg).
addParameter(p, 'cMass', 0, @(x) (x >= 0) && isnumeric(x) && isscalar(x));

% Input mass of reactive tracer (mg).
addParameter(p, 'rMass', 0, @(x) (x >= 0) && isnumeric(x) && isscalar(x));

% Fraction mass recovery of conservative tracer, necessary for proper
% normalization of reactive BTC. If Q is not supplied or if this parameter
addParameter(p, 'fmc', 1, @(x) (x >= 0) && isnumeric(x) && isscalar(x));

% Validate all parameters; extract optional parameters
parse(p, t, c, mymeas, varargin{:});

Q = p.Results.Q;
cMass = p.Results.cMass;

I_BTC =  trapz(t,c); % calculated integral [mg-s/L]

if Q == -inf %If Q not supplied, calculate via dilution gauging.
    if strcmp(mymeas, 'c') && cMass > 0
        Q = cMass / I_BTC; % L/s
        fmc = 1; %Cons. mass recovery assumed 100% with dilution gauging.
    elseif strcmp(mymeas, 'c') && cMass <= 0
        error('Conservative Injection mass -- argument ''cMass'' -- must be > 0 for conservative fit.')

    elseif strcmp(mymeas, 'r')
        error('Discharge must be input to properly normalize reactive BTC.')
    end
elseif Q > 0 % we're ok--do nothing
else
    error ('User-input discharge must be greater than zero.')
end

switch mymeas
    case 'c'
        mc_calc = I_BTC * Q; % recovered conservative mass [mg]
        fmc = mc_calc / cMass; % fraction of conservative mass recovered
        norm_const = cMass * fmc / Q;
        cnorm = c / norm_const; % normalized conservative concentration
        fm = fmc;

    case 'r'
        rMass = p.Results.rMass; % injected reactive mass [mg]
        if rMass <= 0
            error('Reactive tracer mass -- argument ''rMass'' -- greater than zero for normalization.')
        end

        Ir = trapz(t, c);
        mr_calc = Ir * Q;

        if ~exist('fmc', 'var')
            warning('Conservative mass recovery not provided, assuming 100%.')
            fmc = 1;
            norm_const = rMass * fmc / Q;
            cnorm = c / norm_const; % normalized by conservative mass recovery
            fm = mr_calc / rMass * fmc; % fractional conservative mass recovery
        end
end

nout = max(nargout, 1) - 1;
switch nout
    case 1
        varargout{1} = fm;
    case 2
        varargout{1} = fm;
        varargout{2} = Q;
end
