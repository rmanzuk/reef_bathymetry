function [slope,intercept,y_calc,r_squared] = linear_least_squares(x,y)
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
% slope: slope of the linear fit between x and y or 'm' constant in the
% linear equation y = mx + b.
%
% intercept: y-intercept of the linear fit between x and y or 'b' constant
% in the linear equation y = mx + b.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to allow for a possible intercept, we need to add a column of ones to
    % our x vector
    full_x = [ones(length(x),1), x];
    % then calculate
    b = full_x\y;
    slope = b(2);
    intercept = b(1);
    % and give r_squared
    y_calc = full_x*b;
    r_squared = 1 - sum((y - y_calc).^2)/sum((y - mean(y)).^2);
end