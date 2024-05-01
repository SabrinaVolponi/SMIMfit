%===============================================================================
%===============================================================================
%
% This function calculates the Laplace-space memory function, as a function of
% 'sVec' (the Laplace-space variable), for the forward SMIM model.
%
%===============================================================================
%
% Inputs:
%   sVec    := vector of Laplace variables (NOTE: these are complex)
%   options := cell array, containing:
%              (fields used in this function are listed first, as some are
%               simply passing through):
%
%       1) [lambda, tau, W, k_surfrxn, k_subrxn] :=
%                                   Sorption parameters for the memory function
%
%       2) string, indicating the "sticking time" pdf function
%          (for passing to LapPsiFuncTPL)
%          (e.g., 'ADE', 'TPL')
%       3) Standard parameters for the memory function
%          (for passing to LapPsiFuncTPL)
%          (these are [beta, T1, T2], in the case of TPL)
%
% Outputs:
%   memFuncSol := Laplace-space memory function, as a function of 'sVec'
%                 (the Laplace-space variable)
%
%===============================================================================
%===============================================================================

function memFuncSol = memFuncTPL(sVec, options)

if length(options) >= 2  % Sorption/desorption model
    lambda   = options{1}(1);
    tau      = options{1}(2);
    W        = options{1}(3);

    if length(options{1}) > 3 % Added to include tempering KRR20171023
        k_rxn     = options{1}(4 : 5); % first value is k_surfrxn, second value is k_subrxn
        k_surfrxn = k_rxn(1);
        k_subrxn  = k_rxn(2);
    else
        k_surfrxn = 0;
        k_subrxn  = 0;
    end

    lapPsi = LapPsiFuncTPL(sVec + k_subrxn, options{3}, options(2));
    phi = W * lapPsi + (1 - W) ./ (sVec + k_subrxn) / tau;
    % sVec1 is the input Laplace variable for LapPsiFuncTPL
    sVec1 = sVec + lambda .* (1 - phi);
    % NOTE: this may not always be 1--verify this
    tchar = 1; 
    psi = LapPsiFuncTPL(sVec1 + k_surfrxn, [], {'ADE'});
    memFuncSol = tchar * psi ./ (1 - psi) .* sVec;
else
    memFuncSol = 1;
end

end
