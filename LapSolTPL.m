%===============================================================================
%===============================================================================
%
% This function calculates the Laplace-domain solution, as a function of 'sVec'
% (the Laplace-space variable), for the forward SMIM model.
%
%===============================================================================
%
% Inputs:
%   sVec        := vector of Laplace variables (NOTE: these are complex)
%   transParams := vector containing:
%       1) vNorm := normalized velocity
%                   (vNorm := v / L and assumed to be constant)
%       2) DNorm := normalized dispersion coefficient
%                   (Dnorm := D / L^2 and assumed to be constant)
%   options     := cell array, containing:
%              (fields used in this function are listed first, as some are
%               simply passing through):
%
%       1) outletBCtype := Type of outlet (downstream) boundary condition
%                          (NOTE: currently only 'none' is included)
%       2) inletBCfunc  := Function governing the injection at the inlet
%                          (e.g., 'pulse' or specified functional form)
%       3) xSample      := Spatial location where solute is sampled
%                          NOTE: this is defined as a proportion (between 0 and 1)
%                          of the total length (L) of the reach
%       4) inletBCtype  := Type of inlet (upstream) boundary condition
%                          (NOTE: currently only 'none' is included)
%
%       5) Sorption parameters for the memory function
%          (for passing to memFuncTPL)
%          (these are [lambda, tau, W, k_surfrxn, k_subrxn], in the case of TPL)
%       6) String, indicating the "sticking time" pdf function
%          (for passing to LapPsiFuncTPL)
%          (e.g., 'ADE', 'TPL')
%       7) Standard parameters for the memory function
%          (for passing to LapPsiFuncTPL)
%          (these are [beta, T1, T2], in the case of TPL)
%
% Outputs:
%   lapSol := Laplace-space solution to the SMIM forward model
%
%===============================================================================
%===============================================================================

function lapSol = LapSolTPL(sVec, transParams, options)

vNorm        = transParams(1);
DNorm        = transParams(2);
outletBCtype = options{1};
inletBCfunc  = options{2};
xSample      = options{3};
inletBCtype  = options{4};

if vNorm < 0
  error('ERROR: vNorm = %6.5e is negative--Stopping Run', vNorm)
end

if DNorm < 0
  error('ERROR: DNorm = %6.5e is negative--Stopping Run', DNorm);
end

% NOTE: options must be accessed with parentheses, rather than braces, for its
%       output to be a cell array, too
memFunc = memFuncTPL(sVec, options(5 : 7));

switch inletBCfunc
    case {'pulse'}
        % pulse injection/Dirac-delta injection
        LapFuncInletBC = 1;
    otherwise
        % user-defined injection (inletBCfunc is a string providing the
        % functional form)
        LapFuncInletBC = eval(inletBCfunc);
end

alpha = DNorm / vNorm;
expTerm = exp((1 + (4 .* sVec .* alpha) ./ (memFunc * vNorm)).^(1/2) / alpha);

% Note that all results are normalized by their long-time limit in the
% case of a step input, i.e., C/C_inf, where C_inf is the concentration
% for a constant concentration or flux

% NOTE: only the none-none boundary-type case is currently implemented

switch inletBCtype
    case {'none'}
        switch outletBCtype
            case {'none'}	% Robin INLET, Neumann OUTLET
                const = vNorm^2 + 4 * sVec .* DNorm ./ memFunc;
                num = 1 ./ sqrt(const) .*...
                      exp((vNorm - sqrt(const)) / 2 / DNorm * xSample);
                % NOTE NOT SCALED THE SAME!;
                denom1 = 0;
                denom2 = 1;
            otherwise
                error('Inlet and outlet BCs must both be "none"');
        end
end

if strcmp(outletBCtype, 'none')
    denom = 1;
else
    if isinf(abs(expTerm))
        denom = denom1;
    else
        denom = expTerm * denom1 + denom2;
    end
end

lapSol = LapFuncInletBC .* num ./ denom;
