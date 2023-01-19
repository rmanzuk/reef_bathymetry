function [pixel_vals] = track_pixels(bin_inds, photon_ids, track_utmx, track_utmy, satellite_image, grid_x, grid_y, neighborhood_size)
% This function applies a takes the locations where water depth has been
% calculated for an ICESat track and their associated locations and returns
% the pixel values from the associated satellite image.
%
%
% IN: 
%
% bin_inds: discritization indices for all photons, putting them into
% bins based upon distance along the track. From something like:
% along_bins = discretize(gdats.along{track_ind},bin_edges);
%
% photon_ids: column vector of the same length as input_heights with
% indices designating all classes. Output from cluster_photons.m
%
% track_utmx: double vector with the easting coordinates in UTM for photon
% returns.
%
% track_utmy: double vector with the northing coordinates in UTM for photon
% returns.
%
% satellite_image: 3d matrix with the multichannel satellite image from
% which you will be calibrating Beer's Law.
%
% grid_x: 2d mesh grid of utmx coordinates for the satellite_image
%
% grid_y: 2d mesh grid of utmy coordinates for the satellite_image
%
% neighborhood_size: radius of the circular mask used to extract the pixel
% values from the satellite image at each position (for averaging)
%
% OUT: 
%
% pixel_vales: 2d array of the same number of rows as the water depth data
% for the track with associated pixel values for the neighborhood around
% each depth meaurement. Each column represents a band from the satellite
% image.
%
% Written by R. A. Manzuk
% Tuesday, January 17, 2023 at 11:42:52 AM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % we need a pixel size for the neighborhood thing later
    pix_scale = range(grid_x(:))/size(grid_x,2);
    
    % get the mean utmx and utmy position of each bin. we will only use
    % land or surface photon returns be cause bottom photon positions are a
    % little less certain due to refraction. Get rid of noise too.
    [utmx_means, utmy_means] = splitapply(@loc_mean, track_utmx, track_utmy, photon_ids, bin_inds);

    % and find the nearest indices on the utm grid of the satellite image
    pos_utmx_inds = knnsearch(grid_x(1,:)',utmx_means);
    pos_utmy_inds = knnsearch(grid_y(:,1), utmy_means);

    % because we are doing the circular neighborhood thing, probably just
    % easiest to do as a loop for now
    pixel_vals = zeros(numel(pos_utmy_inds), size(satellite_image,3));
    for i = 1:numel(pos_utmy_inds)
        % make a circular mask
        mask = (grid_x - grid_x(1,pos_utmx_inds(i))).^2 + (grid_y - grid_y(pos_utmy_inds(i),1)).^2 <= (neighborhood_size*pix_scale).^2;
        % grab the pixels and reshape into an array with the number of
        % columns matching the number of channels
        these_pixels = reshape(satellite_image(repmat(mask,1,1,size(satellite_image,3))),[],size(satellite_image,3));

        % collapse via the mean and add to the output array
        pixel_vals(i,:) = mean(these_pixels);
    end
end
% function to get the group means of only surface or land pixels

function [mean_utmx,mean_utmy] = loc_mean(utm_x,utm_y,phot_inds)
    % get rid of bottom or noise returns
    good_phots = ismember(phot_inds, [1,3]);
    
    % and calculate the means
    mean_utmx = mean(utm_x(good_phots));
    mean_utmy = mean(utm_y(good_phots));
end