function [adjusted_along] = recalc_along(along, corrected_utmx, corrected_utmy)
% This function just takes the corrected utm positions from refraction
% correction and gives new distances along the track
%
% IN: 
%
% along: vector with the original, uncorrected distances along the track
% for all photons
% 
% corrected_utmx: vector with refraction-corrected easting coordinates for
% all photons.
%
% corrected_utmy: vector with refraction-corrected northing coordinates for
% all photons.
%
% OUT: 
%
% adjusted_along: double vector with new distances along the track for all
% photon returns in meters. The first point is set to a distance of 0.
%
% Written by R. A. Manzuk
% Thursday, January 12, 2023 at 1:58:07 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % use original along data to give a starting point
    [~, starting_ind] = min(along);
    
    % and just pythagorean theorem every point to get a new dist along
    adjusted_along = sqrt((corrected_utmx - corrected_utmx(starting_ind)).^2 +...
        (corrected_utmy - corrected_utmy(starting_ind)).^2);
end