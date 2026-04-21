function beta = qreg_lp(X, y, tau)
% QREG_LP Quantile regression via linear programming.
%
%   beta = qreg_lp(X, y, tau) fits a linear quantile regression model at
%   quantile level tau (0 < tau < 1), returning the coefficient vector
%   beta with an intercept in the first position.
%
%   Inputs:
%       X   - n-by-p predictor matrix (no intercept column; one is added)
%       y   - n-by-1 response vector
%       tau - scalar quantile level in the open interval (0, 1)
%
%   Output:
%       beta - (p+1)-by-1 coefficient vector; beta(1) is the intercept
%
%   Requires the Optimization Toolbox (uses linprog).

    % Add intercept column
    X = [ones(size(X, 1), 1), X];
    [n, p] = size(X);

    % Objective: minimise asymmetric absolute residuals
    f = [zeros(p, 1); tau * ones(n, 1); (1 - tau) * ones(n, 1)];

    % Equality constraints: y = X*beta + u - v, with u, v >= 0
    Aeq = [X, eye(n), -eye(n)];
    beq = y;

    % Bounds: beta unconstrained, residual parts non-negative
    lb = [-inf(p, 1); zeros(2*n, 1)];

    % Solve
    opts = optimoptions('linprog', 'Display', 'none');
    z = linprog(f, [], [], Aeq, beq, lb, [], opts);

    beta = z(1:p);
end
