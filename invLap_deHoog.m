%===============================================================================
%===============================================================================
%
% This function implements the "de Hoog" algorithm for numerically
% calculating an inverse Laplace transform. See:
% https://doi.org/10.1137/0903022
% for the original paper (the equation references below correspond to this
% manuscript).
%
%===============================================================================
%
% Inputs:
%   F    := function handle for Laplace-space function to be transformed
%   tVec := vector of times at which the inverse Laplace transformed
%           function, f(t) will be evaluated
%
% Outputs:
%   invLap := inverse Laplace-transformed function, f(t)
%
%===============================================================================
%===============================================================================

function invLap = invLap_deHoog(F, tVec)

% make sure tVec is a column
if ~iscolumn(tVec)
    tVec = tVec';
end

[tol, M] = getTolerance('invLap_deHoog');

%  alpha is the real part of the rightmost pole or singularity
%     mjs(FIXME?): this appears to be commonly fixed to be zero (in CTRW
%         toolbox and the CPP code referenced below.
%     mjs(FIXME?): should this be moved to getTolerance()?
alpha = 0;
% number of terms in Taylor approximation
ns = 2 * M + 1;

% group the time vector by orders of magnitude and approximate separately
%     on each grouping--this improves the approximation.
%     see Kuhlman, "Review of inverse Laplace transform
%     algorithms for Laplace-space numerical approaches," Numerical
%     Algorithms, 2013, for a discussion of this principle.
mags = floor(log10(tVec));
minMag = min(mags);
numMags = max(mags) - minMag + 1;
tCell = cell(numMags, 1);
for i = 1 : numMags
    tCell{i} = tVec(mags == minMag + i - 1);
end

% this stores the individual solutions for each time grouping
invLapCell = cell(numMags, 1);

for i = 1 : numMags

    tVec = tCell{i};
    nt = length(tVec);

%     scaled time--2 appears to be the typical choice
    T = 2 * max(tVec);
%     note: natural log seems to be the appropriate choice here, based on
%     testing. See also CPP and Python versions of this algorithm that
%     employ natural log:
%         https://www.civil.uwaterloo.ca/jrcraig/pdf/LaplaceInversion.cpp
%         http://mpmath.org/doc/current/calculus/inverselaplace.html
    gamma = alpha - log(tol) / (2 * T);

%     vector of Laplace variables (e.g., 's')
    sVec = gamma + (1i * pi * (0 : 2 * M)') / T;

%     evaluate F(s) at the points defined above
    Fs = F(sVec);

%     pre-allocate some vectors
    e = zeros(ns, M + 1);
    q = zeros(ns - 1, M);
    d = zeros(ns, 1);
    A = complex(zeros(nt, ns + 1));
    B = complex(ones(nt, ns + 1));

%     Eq. (20) in de Hoog--dividing the first entry of F(s) by 2 is
%     referenced below Eq. (21)
    Fs(1) = Fs(1) / 2;
    q(:, 2) = Fs(2 : ns) ./ Fs(1 : ns - 1);
    for r = 2 : M + 1
        mr = 2 * (M - r + 1) + 1;
        e(1 : mr, r) = q(2 : mr + 1, r) -...
                       q(1 : mr    , r) +...
                       e(2 : mr + 1, r - 1);
        if r < M + 1
            rq = r + 1;
            mr = 2 * (M - rq + 1) + 3;
            q(1 : mr - 1, rq) = q(2 : mr    , rq - 1) .*...
                                e(2 : mr    , rq - 1) ./...
                                e(1 : mr - 1, rq - 1);
        end
    end

%     continued fraction coefficients--unnumbered Eq. on p. 362 of de Hoog
    d(1)              = Fs(1);
    d(2 : 2 : ns - 1) = -q(1, 2 : M + 1);
    d(3 : 2 : ns)     = -e(1, 2 : M + 1);

%     initialize A for recurrence and define z (text following Eq. (21) in de Hoog)
%     note that B is already initialized above
    A(:, 2) = d(1);
    z = exp(1i * pi * tVec / T);

%     Eq. (21) in de Hoog
    for n = 3 : ns
        A(:, n) = A(:, n - 1) + d(n - 1) .* z .* A(:, n - 2);
        B(:, n) = B(:, n - 1) + d(n - 1) .* z .* B(:, n - 2);
    end

%     remainder terms--Eq (23) in de Hoog
    h2M  = (1 + (d(ns - 1) - d(ns)) .* z) / 2;
    R2M   = h2M .* ((1 + d(ns) .* z ./ h2M.^2).^0.5 - 1);

%     Eq. (24) in de Hoog
    A(:, ns + 1) = A(:, ns) + R2M .* A(:, ns - 1);
    B(:, ns + 1) = B(:, ns) + R2M .* B(:, ns - 1);

%     Eq. (25) in de Hoog--f(t) for this time grouping
    invLapCell{i} = (exp(gamma * tVec) / T) .* real(A(:, ns + 1) ./ B(:, ns + 1));

end

% full f(t) solution for all times in tVec
invLap = cell2mat(invLapCell);

end
