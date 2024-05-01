%===============================================================================
%===============================================================================
%
% This function creates the objective function to be used in the minimization
% problem that fits the SMIM forward model to data.
% This function takes inspiration from: Kelly, James F., et al.
% "FracFit: A robust parameter estimation tool for fractional calculus models."
% Water Resources Research 53.3 (2017): 2559-2567.
%
%===============================================================================
%
% Inputs:
%   params := parameters sent to the pdf function
%   K_mass := currently not used, always set to 1
%   data   := struct file containing fields for
%       - observation times (.tobs)
%       - observed concentrations (.cobs)
%       - times from conservative fit, if completed (.tcfit)
%       - concentrations from conservative fit, if completed (.ccfit)
%       - length of domain (.L)
%       - duration of injection at inlet boundary (.injectDuration)
%   pdf_function := function handle for the SMIM forward model,
%                   set in SMIMfit()
%
% Outputs:
%   f := objective function to be minimized
%        (weighted absolute error between data and forward model)
%
%===============================================================================
%===============================================================================

function f = objFunctionSMIM(params, K_mass, data, pdf_function)

cobs = data.cobs;
N_samples = length(cobs);
c_fit = K_mass .* pdf_function(params, data);
if(sum(isnan(c_fit)) ~= 0 )
    f = repelem(1e16, numel(cobs));
    return
else

% KRR 20190920: second scalefactor calculates the expected model recovery, which corrects
% for mass that will have passed in a window outside tobs. This correction
% is currently necessary because we are assuming complete mass recovery in
% the interval tobs.
if isfield(data, 'tcfit') && isfield(data, 'ccfit')
    scalefactor = trapz(data.tcfit, data.ccfit)
    %only populated when performing reactive fits
else
    % else scalefactor 2 will be adjusting during conservative fit
    scalefactor = trapz(data.tobs, c_fit);

end

c_fit = c_fit / scalefactor;

% weights are from Kelly et al., "FracFit" paper
wgts = 1.0 ./ sqrt(N_samples * K_mass .* cobs);

% this is the current algorithm for handling zeros--the old one used interpolation
wgts(isinf(wgts)) = 0;


%Calculate weighted mean absolute error.
f = wgts .* abs(cobs - c_fit);
end
end

