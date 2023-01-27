function [r_squareds, best_x] = variance_reduction(x_observations, y_observations)
% This function performs a sequential regression through a set of x
% variables that may predict a y variable. At each step it identifies the
% variable that most greatly adds to the predictive power of the regression
% (through r^2), creating a sequential list of the variables and the
% r_squareds of all regressions.
%
% IN: 
% 
% x_observations: matrix where each column is a different x variable's
% observations
%
% y_observations: single column matrix with the y observations to be
% predicted.
%
% OUT: 
%
% r_squareds: square 2D matrix with the same number of rows and columns as
% the number of variables (columns) in the input x_observations. Each
% column is the r_squared associated with each variable for that particular
% round of regression. 1st column is the start, so it's just the r_squared
% for each variable on its own. 0s indicate the variable was already added
% to the existing regression and is no longer being assessed for additive
% power.
%
% best_x: sequence of variable indices that indicate the order in which
% variables are added to the regression based upon maximum contribution to
% r^2
%
% Written by R. A. Manzuk
% Friday, January 27, 2023 at 11:12:01 AM
%
% Edited: Friday, January 27, 2023 at 11:12:06 AM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % we're going to return all of the r_squared values that we calculate
    % iteratively here. So we'll start with an empty matrix and fill it up.
    r_squareds = zeros(size(x_observations,2));

    % and a matrix to fill with the best x variable added at each iteration
    best_x = zeros(1, size(x_observations,2));
    for i = 1:size(x_observations,2)
        % start by deciding which variables haven't been attributed to the
        % regression yet
        not_regressed = setdiff(1:20,best_x);
        
        % and make a matrix of x observations that have already been
        % determined to be part of the best fit
        prev_x = x_observations(:,best_x(best_x ~= 0));
        for j = not_regressed
            [~, ~, ~, r_squareds(j,i)] = linear_least_squares([prev_x, x_observations(:,j)], y_observations);
        end
        [~, best_x(i)] = max(r_squareds(:,i));
    end
end