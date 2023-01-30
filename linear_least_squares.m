function [coefficients,y_calc,r_squared] = linear_least_squares(x,y)
% Matlab has a lot of regression functions and I don't really like them. So
% good to have my own simple linear least squares ready to go when needed.
%
% INPUT 
%
% x: vector with observations for the variable to be considered as 'x' in
% the linear equation y = mx + b.
%
% y: vector with observations for the response or 'y' variable in the
% linear equation y = mx + b.
%
% OUTPUT
%
% coefficients: coeffiecients  of the linear fit between x and y. The first
% coefficient is the 'intercept' and the following ones correspond to each
% variable from the input x observations.
%
% y_calc: vector with calculated y values using the slope, intercept, and
% input x vector.
%
% r_squared: proportion of the y variable that is explained by the x
% variable. Typical 'goodness of fit' metric. Calculated with the formula:
% R^2 = 1 - (sum of suqred residuals / total sum of squares)
%
% 
% Written by R.A. Manzuk
% 11/28/2022
% Edited: Monday, January 30, 2023 at 4:09:15 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to allow for a possible intercept, we need to add a column of ones to
    % our x vector
    full_x = [ones(length(x),1), x];
    % then calculate
    coefficients = full_x\y;
    % and give r_squared
    y_calc = full_x*coefficients;
    r_squared = 1 - sum((y - y_calc).^2)/sum((y - mean(y)).^2);
end