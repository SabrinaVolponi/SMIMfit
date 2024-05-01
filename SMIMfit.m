%===============================================================================
%===============================================================================
%
% This function is the high-level script that conducts the fit of the forward
% SMIM model to data.
%
%===============================================================================
%
% Inputs:
%   mymeas   := measurement (conservative or reactive)
%   t        := time vector corresponding to concentration data
%   c        := concentration vector (mg/L)
%   varargin := optional inputs (must be specified as " 'name', value " pairs):
%       - model_type     := Which model you're using. Only programmed for truncated
%                           power law at the moment (default TPL).
%       - L              := Length of the reach where measures are being made
%                           (default 50 m).
%       - injectDuration := duration of the injection. Exclude this name/value pair
%                           or set to zero if injection is a pulse.
%       - t_end          := latest time you wish to use in the model fit. Set this
%                           value when tailing behavior is unphysical due to,
%                           e.g., return to background conditions.
%       - M              := Output of a previous run (i.e., the struct output).
%                           Use to quickly specify guesses and limits of
%                           parameter fits (if non-default) or to supply modeled
%                           time series of conservative tracer when fitting a
%                           reactive time series.
%
% Outputs:
%   M := struct containing results of the fit, containing:
%       - params_fit      := fitted parameters (i.e., the results)
%       - opts            := fit options, including the parameter guesses and
%                            upper/lower bounds
%       - obj_fcn_result  := objective function residual vector from best fit
%       - SSE             := sum of squared errors from best fit
%       - fitdate         := time and date of final fit
%       - allfits         := parameter array containing fitted parameters for each fit
%       - obj_fcn_allfits := norm of objective function residual for each fit
%       - algout          := struct containing output from the optimization function
%       - cXfit           := best fit concentration vector, where 'X' is 'c' or
%                            'r', depending on whether this is for a conservative
%                            or reactive fit
%       - tXfit           := time vector corresponding to cXfit, where 'X' is
%                            'c' or 'r', depending on whether this is for a
%                            conservative or reactive fit
%       - kfits           := only generated when performing reactive fit,
%                            contains structs for the type of reactive fit that
%                            is conducted (i.e., kboth, ksurf, ksub), and these
%                            structs will contain the above fields for the reactive fit
%
%===============================================================================
%===============================================================================

function M = SMIMfit(t, c, mymeas, varargin)

% 2019/11/27 - Removed newguess feature and cleaned up parameter_guess
% inputs. Also, updated CheckGuess to properly report errors.
% 2019/12/17 - Added optional parameter 'injectDuration' to allow for step
% injections of a specified duration.
% 2020/02/02 - Add field L to struct 'data' to allow for different reach sizes


%% Function Setup
p = inputParser;

% Error checking functions
check_t = @(x) validateattributes(x, {'double'}, {'increasing', 'nonempty', 'nonnegative'});
check_c = @(x)  validateattributes(x, {'double'}, {'nonempty', 'nonnegative'});
checkMeas = @(x) any(validatestring(x, {'c', 'r'}));
checkModel = @(x) any(validatestring(x, {'TPL'}));


% Required input parameters
addRequired(p, 't', check_t); % time vector
addRequired(p, 'c', check_c); % concentration vector
addRequired(p, 'mymeas', checkMeas); % measurement, must be 'c' (conservative) or 'r' (reactive)

% Optional input parameters; variables are set to defaults unless called as
% a 'name',value argument pair in the function input
addParameter(p, 'model_type', 'TPL', checkModel); % Which model you're using. Only programmed for truncated power law at the moment (default TPL).
addParameter(p, 'L', 50, @(x) (x > 0) && isnumeric(x) && isscalar(x)); % Length of the reach where measures are being made (default 50 m).
addParameter(p, 'injectDuration', 0, @(x) (x >= 0) && isnumeric(x) && isscalar(x)); % duration of the injection. Exclude this name/value pair or set to zero if injection is a pulse.

% latest time you wish to use in the model fit. Set this value when tailing behavior
% is unphysical due to, e.g., return to background conditions.
addParameter(p, 't_end', inf, @(x) (x > 0) && isnumeric(x) && isscalar(x));
addParameter(p, 'params_guess', [], @(x) or(length(x) == 6, length(x) == 8));
addParameter(p, 'params_upper', [], @(x) or(length(x) == 6, length(x) == 8));
addParameter(p, 'params_lower', [], @(x) or(length(x) == 6, length(x) == 8));



% Output of a previous model run. Use to quickly specify guesses and limits
% of parameter fits (if non-default) or to supply modeled time series of
% conservative tracer when fitting a reactive time series.
addParameter(p, 'M', struct, @(x) isstruct(x));

% Validate all parameters; extract optional parameters
parse(p, t, c, mymeas, varargin{:});

model_type     = p.Results.model_type;
L              = p.Results.L;
t_end          = p.Results.t_end;
M              = p.Results.M;
injectDuration = p.Results.injectDuration;

if(~isempty(p.Results.params_guess))
    M.params_guess = p.Results.params_guess;
    M.params_upper = p.Results.params_upper;
    M.params_lower = p.Results.params_lower;

end


%% Parse the time Vector to ensure it is valid

tic %start timer

% Remove datapoints later than t_end.
if ~isempty(t_end)
    idx = find(t > t_end);
    t(idx) = [];
    c(idx) = [];
    clear idx
end

if t(1) == 0
    t = t(2 : end);
    c = c(2 : end);
end

% only fit where concentration values are strictly positive
idx = find(c <= 0);
t(idx) = [];
c(idx) = [];

% Find the time of peak C (for initial guess of velocity)
imax = find(c == max(c), 1, 'first');

data.tobs = t;
data.cobs = c;
data.L = L;


%%

% create this text variable for printing to screen
if strcmp(mymeas, 'c')
    mymeas_text           = 'Conservative';
elseif strcmp(mymeas,'r')
    mymeas_text           = 'Reactive';
else
    error('Field ''mymeas'' must be set to ''c'' or ''r''');
end

fprintf('\nStarting %s %s fits for L = %3.2f m\n', mymeas_text, model_type, L)

%% Initial parameter guess

% The options here are:
%     1) Conservative: use the default guess, found in the
%     initialGuessTPL() function, or pass your own guess in via the
%     provided data struct, M
%     2) Reactive: use the initialGuessTPL() function to modify the guess for
%     the reactive  parameters

if strcmp(mymeas, 'c')
    if isfield(M, 'params_guess') % if provided, use the input from the struct M
        params_guess   = M.params_guess;
        params_lower   = M.params_lower;
        params_upper   = M.params_upper;


    else
        switch model_type
            case 'TPL'
                [params_guess, params_upper, params_lower] =...
                       initialGuessTPL(mymeas, 'length', L, 'CmaxTime', t(imax));
        display(params_lower)
        display(params_upper)
        display(params_guess)
        end
    end
%     Check they are valid guesses (see function checkGuess at end of file)
    checkGuess(params_guess, model_type, mymeas)
    checkGuess(params_upper, model_type, mymeas)
    checkGuess(params_lower, model_type, mymeas)
end


if strcmp(mymeas, 'r') % only fit for reaction parameters
    if isfield(M, 'params_guess') % if provided, use the input from the struct M
        params_guess   = M.params_guess;
        params_lower   = M.params_lower;
        params_upper   = M.params_upper;
    else
        switch model_type
            case 'TPL'
                [params_guess, params_upper, params_lower] =...
                                                        initialGuessTPL(mymeas, 'consM', M);

        end
    end
%     Check they are valid guesses (see function checkGuess at end of file)
    checkGuess(params_guess, model_type, mymeas)
    checkGuess(params_upper, model_type, mymeas)
    checkGuess(params_lower, model_type, mymeas)
    
    % check to see if conservative fit is complete
    if isfield(M, 'ccfit') && isfield(M, 'tcfit')
        data.tcfit = M.tcfit;
        data.ccfit = M.ccfit;
    else
        error('Conservative fits must be complete before performing reaction fits')
    end

    % store transport parameters so we don't change them during every reactive
    % fit iteration
    params_transport = params_guess(1 : end - 2);


    % option to fit just the surface just the subsurface, or both
        % e.g., ktypes = {'ksurf', 'ksub', 'kboth'};
    ktypes = {'kboth'};


    % The guesses and upper/lower bounds will likely need to change based
    % on the specific dataset.
    kfits = struct;
    k_guess = params_guess(end - 1 : end);
    k_upper = params_upper(end - 1 : end);
    k_lower = params_lower(end - 1 : end);

    for jj =1 : numel(ktypes)
        kfits.(ktypes{jj}).k_guess = k_guess;
        kfits.(ktypes{jj}).k_upper = k_upper;
        kfits.(ktypes{jj}).k_lower = k_lower;

        switch ktypes{jj}
            case 'ksurf'
                kfits.(ktypes{jj}).k_guess(2) = 0;
                kfits.(ktypes{jj}).k_upper(2) = 0;
                kfits.(ktypes{jj}).k_lower(2) = 0;

            case 'ksub'
                kfits.(ktypes{jj}).k_guess(1) = 0;
                kfits.(ktypes{jj}).k_upper(1) = 0;
                kfits.(ktypes{jj}).k_lower(1) = 0;

            case 'kboth'
        end
    end
end

data.injectDuration = injectDuration;
 %create function handles for objective and pdf functions
[obj_function, pdf_function] = createModel(model_type);

%% Run the fit

switch mymeas
    case 'c'

%         run the fit
        M = runFit(mymeas, obj_function, data, pdf_function,...
                   params_guess, params_upper, params_lower);

%         write some additional info to the data struct
        M.tcfit = t;
        M.cobs = c;
        
%         print a summary to the screen
        fprintf('Fit complete');
        
    case 'r'
        for jj = 1 : numel(ktypes)
            k_guess      = kfits.(ktypes{jj}).k_guess;
            k_upper      = kfits.(ktypes{jj}).k_upper;
            k_lower      = kfits.(ktypes{jj}).k_lower;
            params_guess = k_guess;
            params_upper = k_upper;
            params_lower = k_lower;

%             run the fit
            kfits.(ktypes{jj}) = runFit(mymeas, obj_function, data, pdf_function,...
                                        params_guess, params_upper, params_lower,...
                                        'transParams', params_transport);

%             write some additional info to the data struct
            kfits.(ktypes{jj}).trfit       = t;

%             print a summary to the screen
            fprintf(['\nPrintout of all reactive fit iterations for experiment, ',...
                     'fits for %s.\n'], ktypes{jj});
            fprintf(['Conservative params: U %2.2f m/s, D %2.2e cm/s ',...
                     'Lambda %2.3f, Beta %2.3f, logT %2.2f.\n'],...
                    params_transport(1), params_transport(2), params_transport(3), ...
                    params_transport(4), params_transport(6));

            for n_iter = 1 : size(kfits.(ktypes{jj}).params_fit, 1)
                fprintf('iter %02u: k_surf %3.3e, ksub %3.3e, SSE %3.3e.\n',...
                        n_iter, kfits.(ktypes{jj}).params_fit(n_iter, 1),...
                        kfits.(ktypes{jj}).params_fit(n_iter, 2),...
                        kfits.(ktypes{jj}).obj_fcn_allfits(n_iter))
            end
        end
        M.kfits = kfits;
end

fprintf('\nCompleted %s %s fits. Time = %f min.\n', mymeas_text, model_type, toc / 60)
end

function checkGuess(x, model_type, meas)
lx = length(x);
switch model_type
    case 'TPL'
        if (lx ~= 6 && strcmp(meas, 'c')) || (lx ~=8 && strcmp(meas, 'r'))
            error(['Guess and bounds vectors for TPL model must be length 6 ',...
                   '(conservative) or 8 (reactive).'])
        end
end
end
