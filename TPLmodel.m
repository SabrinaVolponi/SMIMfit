%===============================================================================
%===============================================================================
%
% This function runs a forward SMIM TPL (truncated power law) model, returning a
% concentration vector that is a function of space.
%
%===============================================================================
%
% Inputs:
%   allParams := parameter vector, holding:
%       1) v      := velocity
%       2) D      := diffusion coefficient
%       3) Lambda := average "sticking rate" (a general sticking-time pdf parameter)
%       4) beta   := parameter for TPL sticking-time pdf
%       5) t1     := parameter for TPL sticking-time pdf
%       6) t2     := parameter for TPL sticking-time pdf
%   data      := struct file containing fields for:
%       - observation times (.tobs)
%       - observed concentrations (.cobs)
%       - times from conservative fit, if completed (.tcfit)
%       - concentrations from conservative fit, if completed (.ccfit)
%       - length of domain (.L)
%       - duration of injection at inlet boundary (.injectDuration)
%
% Outputs:
%   c := concentration vector
%
%===============================================================================
%===============================================================================

function c = TPLmodel(allParams, data)

% note: this is the duration of the injection, not the time step length (dt)
injectDuration = data.injectDuration;

if injectDuration > 0
    inletBCfunc = sprintf('(1 - exp(-%f * u)) / u', injectDuration);
else
    inletBCfunc = 'pulse';
end

v = allParams(1);
D = allParams(2);
Lambda = allParams(3);
paramsTPL = allParams(4 : 6);

switch length(allParams)
    case 6
        k_surfrxn = 0;
        k_subrxn = 0;
    case 8
        k_surfrxn = allParams(7);
        k_subrxn = allParams(8);
    otherwise
        error('TPLmodel input ''allParams'' must have length 6 (conservative) or 8 (reactive)')
end

t = data.tobs;

L = data.L;
vnorm = v / L;
Dnorm = D / L^2;

transParams = [vnorm Dnorm];
options = {'none', inletBCfunc, 1, 'none', [Lambda 1 1 k_surfrxn k_subrxn],...
           'TPL', paramsTPL};

fun = @(s) LapSolTPL(s, transParams, options);

c = invLap_deHoog(fun, t);


if numel(allParams) == 8 && sum(allParams(end - 1 : end)) > 0
    % Use this normalization for the reactive case, i.e., when reaction
    % variables are declared and these variables are greater than zero. This
    % normalization is equal to the expected mass recovered in the
    % observation time interval when fitting by the CTRW with conservative
    % transport.
    
    scalefactor = trapz(data.tcfit, data.ccfit);

    if injectDuration > 0 % Normalized by total mass!
        c = c / scalefactor * vnorm / injectDuration;
    else
        c = c / scalefactor * vnorm;
    end
else
    if injectDuration > 0
        c = c * vnorm / injectDuration;
    else
        c = c * vnorm;
    end
end

end
