%===============================================================================
%===============================================================================
%
% This function provides parameters (related to approximation accuracy or
% tolerance) for numerical methods.
%   NOTE: currently only invLap_deHoog() and GaussLaguerre() call this function
%
%===============================================================================
%
% Inputs:
%   fName := name of the calling function--determines what parameters will be provided
%
% Outputs:
%   varargout := holds the relevant parameters for the calling function:
%       - invLap_deHoog:
%           1) tol := the specified accuracy of the method
%           2) M   := degree of Taylor expansion to be used
%       - GaussLaguerre:
%           1) polyDegree := degree of Laguerre polynomial to be used for
%                            Gauss-Laguerre quadrature
%           2) alpha      := generalized Laguerre polynomial parameter--
%                            alpha = 0 corresponds to simple Laguerre polynomials

%===============================================================================
%===============================================================================

function [varargout] = getTolerance(fName)

p = inputParser;

% error checking function
checkfName = @(x) any(validatestring(x, {'invLap_deHoog', 'GaussLaguerre'}));

% calling function name
addRequired(p, 'fName', checkfName);

parse(p, fName);

switch fName
    case{'invLap_deHoog'}
%         specified accuracy of the method
        tol = 1e-10;
%         degree of Taylor expansion
        M   = 25;

        varargout{1} = tol;
        varargout{2} = M;
    case{'GaussLaguerre'}
%         degree of Laguerre polynomial to be used for quadrature
        polyDegree = 8;
%         generalized Laguerre polynomial parameter--alpha = 0 is simple version
        alpha = 0;

        varargout{1} = polyDegree;
        varargout{2} = alpha;
end

end
