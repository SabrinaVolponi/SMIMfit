%===============================================================================
%===============================================================================
%
% This function creates the forward SMIM model to be used by specifying the
% pdf/cdf function and the objective function to be minimized when fitting the
% forward model to data in SMIMfit().
%
%===============================================================================
%
% Inputs:
%   modelName := string specifying the forward model to be employed
%       NOTE: TPL (truncated power law) is currently the only model available
%
% Outputs:
%   obj_function := a function handle for the objective function to be minimized
%   xdf_function := a function handle for the pdf or cdf function to be employed
%                   as the forward model
%
%===============================================================================
%===============================================================================

function [obj_function, xdf_function] = createModel(modelName)

switch modelName
    case 'TPL'
        obj_function = @objFunctionSMIM;
        xdf_function = @TPLmodel;
end

end
