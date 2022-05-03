% Thanks to Henry for providing this code
% Determines the LCP error. This measure of error is unit consistent, i.e. 
% the unit is in J. We are measuring the energy error.
%
% The system has the following dimensions:
%   n = number of bodies (parts) 
%   k = number of constraints 
%   m = number of contraint rows 
%
% Inputs:
%   w   - m-by-1 vector that represents the constrained velocity, where
%         A*x - b = w = w_+ w_-. w_+ and w_- is determined by sign of w
%
%   x   - m-by-1 vector of contraint forces (lambda)
%
%	A_diag_vec  - the m-by-1 vector that consists of diagonal entris of A.
%   lower       - m-by-1 lower bound vector, -inf for unbounded constraints
%   upper       - m-by-1 upper bound vector, +inf for unbounded constraints
function error = LCP_error( w, x, A_diag_vec, lower, upper )
    % Make sure the dimension matches
    assert(length(x) == numel(w) && numel(x) == numel(lower) && numel(x) == numel(upper));

    % Compute the violation of the upper and lower bound
    delta_x_u = max(x - upper, 0);
    delta_x_l = max(lower - x, 0);
    
    % Compute feasible component x0
    x0 = max( lower, min( x, upper ) );
    
    % Compute upper impulse energy error
    delta_e_x_u = 0.5 * A_diag_vec .* (delta_x_u .* delta_x_u);
    
    % Compute lower impulse energy error
    delta_e_x_l = 0.5 * A_diag_vec .* (delta_x_l .* delta_x_l);
    
    % Compute w_+ and w_-
    % Note: For component i, if w(i) > 0, w_plus(i) = w(i), w_minus(i) = 0.
    % Likewise if w(i) < 0, w_minus(i) = -w(i), w_plus(i) = 0
    w_plus = max(w, 0);
    w_minus = -min(w, 0);
    
    % Compute the saturation of the upper and lower bound
    sigma_l = (x0 + delta_x_u) - lower;
    sigma_u = upper - (x0 - delta_x_l);
    
    % Compute the positive velocity energy error
    delta_e_w_plus = min( ((1 ./ (2 * (A_diag_vec + 1e-10))) .* (w_plus .* w_plus)), ...
                          (0.5 * A_diag_vec .* (sigma_l .* sigma_l)) );
    
    % Compute the negative velocity energy error                  
    delta_e_w_minus = min( ((1 ./ (2 * (A_diag_vec + 1e-10))) .* (w_minus .* w_minus)), ...
                          (0.5 * A_diag_vec .* (sigma_u .* sigma_u)) );
                      
    % Compute the energy error per constraint
    delta_e = max( max( abs(delta_e_x_u), abs(delta_e_x_l) ), ...
                   max( abs(delta_e_w_plus), abs(delta_e_w_minus) ) );
    
    % Compute the LCP error (the total error of the system)
    error = norm(delta_e, 1);
end