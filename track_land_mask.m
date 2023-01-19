function [track_mask] = track_land_mask(track_utmx, track_utmy, land_mask, grid_x, grid_y)
% This function applies a land image mask (true for ocean, false for land)
% to icesat photon return track to give a logical mask for the track to get
% rid of land
%
%
% IN: 
%
% track_utmx: double vector with the easting coordinates in UTM for photon
% returns.
%
% track_utmy: double vector with the northing coordinates in UTM for photon
% returns.
%
% land_mask: 2d logical where false is land and true is sea in the
% corresponding satellite imagery.
%
% grid_x: 2d mesh grid of utmx coordinates for the mask
%
% grid_y: 2d mesh grid of utmy coordinates for the mask
%
% OUT: 
%
% track_mask: logical vector of the same length as the input track utmx
% coordinates where false represents land photon returns to be masked out
%
% Written by R. A. Manzuk
% Friday, January 13, 2023 at 8:19:02 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % first thing is to profile the mask, and utm mesh grids based upon the
    % track
    utmx_start = track_utmx(1);
    utmx_end = track_utmx(end);
    utmy_start = track_utmy(1);
    utmy_end = track_utmy(end);

    % and use those starting locs to find the closest pixels in the utm grids
    [~,utmx_start_closest] = min(abs(grid_x(1,:) - utmx_start));
    [~,utmx_end_closest] = min(abs(grid_x(1,:) - utmx_end));
    [~,utmy_start_closest] = min(abs(grid_y(:,1) - utmy_start));
    [~,utmy_end_closest] = min(abs(grid_y(:,1) - utmy_end));
    
    % now take the profiles
    mask_profile = improfile(land_mask,[utmx_start_closest,utmx_end_closest],[utmy_start_closest,utmy_end_closest]);
    % and the associated locations with the profile
    utmx_profile = improfile(grid_x,[utmx_start_closest,utmx_end_closest],[utmy_start_closest,utmy_end_closest]);
    utmy_profile = improfile(grid_y,[utmx_start_closest,utmx_end_closest],[utmy_start_closest,utmy_end_closest]);

    % now we need to do a nearest neighbor search between the track utm
    % coordinates and those from the pixel grid to assign every photon
    % return a position on the mask profile
    track_mask_inds = knnsearch([utmx_profile,utmy_profile],[track_utmx,track_utmy]);
    
    % boom, just have to nearest neighbors to index the mask profile.
    track_mask = logical(mask_profile(track_mask_inds));
end