%===============================================================================
%===============================================================================
%
% This function calculates the Laplace-space 'psi,' or "sticking time,"
% function, as a function of 'sVec' (the Laplace-space variable), for the
% forward SMIM model.
%
%===============================================================================
%
% Inputs:
%   sVec    := vector of Laplace variables (NOTE: these are complex)
%   param   := vector, containing standard parameters for the memory function:
%              [beta, T1, T2]
%   pdfFunc := string, indicating the "sticking time" pdf function (for passing to LapPsiFuncTPL)
%              (e.g., 'ADE', 'TPL')
%
% Outputs:
%   LapPsiSol := Laplace-space 'psi,' or "sticking time," function, as a
%                function of 'sVec'(the Laplace-space variable)
%
%===============================================================================
%===============================================================================

function LapPsiSol = LapPsiFuncTPL(sVec, param, pdfFunc)

switch pdfFunc{1}
    case {'ADE'}
        tchar = 1; % the exponential psi of the ADE model is not used directly
        LapPsiSol  = 1 ./ (1 + tchar .* sVec);
    case {'TPL'}
        beta    = param(1);
        T1      = 10^param(2);% NOTE: for the TPL model, log_10(T1) and log_10(T1) are
        T2      = 10^param(3);%       passed, hence the exponentiation
        invTau2 = T1 ./ T2;
        x       =  invTau2 + (T1 .* sVec);  % for small invTau2, x <= 1 when T1 > u, or for large times relative to T1?
        denom   = (invTau2.^beta .* exp(invTau2) .* igamma(-beta, invTau2)).^(-1);

        % FIXME(mjs): Is this doing what it's supposed to be doing?
        % note that because x here (the Laplace variable vector) is complex,
        % this inequality only applies to the real portion of x
        if any(x > 1)

            idx = x <= 1;
            LapPsiSol = zeros(size(x));

            % this is the portion corresponding to x <= 1
            LapPsiSol(idx) = denom(idx) .* x(idx).^beta .* exp(x(idx)) .* igamma(-beta, x(idx));

            warning(['This case, employing Gauss-Laguerre quadrature, has not ',...
            'been fully tested--please email Michael Schmidt with your data ',...
            'and/or results for testing: mschmi23@nd.edu'])
            % % this is the portion corresponding to x > 1, employing Gauss-Laguerre quadrature
            [polyDegree, alpha] = getTolerance('GaussLaguerre');
            [pts, wts] = GaussLaguerre(polyDegree, alpha);
            xGL = x(~idx);
            intGL = zeros(length(xGL), 8);

            for i = 1 : polyDegree
                fp = (1 + pts(i) ./ xGL).^( -1 - beta);
                intGL(:, i) = wts(i) .* fp;
            end
            sum_intGL = sum(intGL, 2);
            LapPsiSol(~idx) = denom(~idx) ./ xGL .* sum_intGL;

        else % all entries of x <= 1
            LapPsiSol = denom .* x.^beta .* exp(x) .* igamma(-beta, x);
        end
end

end
