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
% values from the satellite image at each position (for averaging). If
% neighborhood size is zero, the function just gets the single pixel value
% at the specified location.
%
% OUT: 
%
% pixel_vales: 2d array of the same number of rows as the water depth data
% for the track with associated pixel values for the neighborhood around
% each depth meaurement. Each column represents a band from the satellite
% image. some places don't have imagery, so we leave zeros there.
%
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
    % do in for loop for now. Not that slow anyways
    utmx_means = zeros(max(bin_inds),1);
    utmy_means = zeros(max(bin_inds),1);
    for i = 1:numel(utmy_means)
        these_inds = bin_inds == i;
        [utmx_means(i), utmy_means(i)] = loc_mean(track_utmx(these_inds), track_utmy(these_inds), photon_ids(these_inds));
    end

    % and find the nearest indices on the utm grid of the satellite image
    pos_utmx_inds = knnsearch(grid_x(1,:)',utmx_means);
    pos_utmy_inds = knnsearch(grid_y(:,1), utmy_means);

    if neighborhood_size == 0
        % in this case, we just grab the pixels
        % zeros array to fill in pixel values where we have imagery
        pixel_vals = zeros(numel(pos_utmy_inds), size(satellite_image,3));

        % check which neighborhoods are complete in the image
        valid_x = pos_utmx_inds <= size(grid_x,2) & pos_utmx_inds > 0;
        valid_y = pos_utmy_inds <= size(grid_x,1) & pos_utmy_inds > 0;
        has_imagery = all([valid_x,valid_y],2);
       
        % just use linear indices for speed and simplicity. As a column vector
        % for now. We can smartly reshape later
        linear_inds = sub2ind(size(grid_x), pos_utmy_inds, pos_utmx_inds);

        % make the image into a 2d matrix
        col_image = reshape(satellite_image, [], size(satellite_image,3));

        % and just index the proper pixels
        pixel_vals(has_imagery,:) = col_image(linear_inds(has_imagery), :);

    else
        % to quickly get the pixel values, we'll have to create a set of pixels
        % in the circle neighborhoods around each pixel. so we'll start by just
        % making one neighborhood that we'll translate to all the other
        % locations. 
        % start with a mask for the first neighborhood
        mask = (grid_x - grid_x(1,pos_utmx_inds(1))).^2 + (grid_y - grid_y(pos_utmy_inds(1),1)).^2 <= (neighborhood_size*pix_scale).^2;
    
        % get the indices of all pixels within that mask
        [y_inds_init,x_inds_init] = find(mask == true);
    
        % define the displacement between all the sampled points and the first
        % point
        x_displacements = pos_utmx_inds - pos_utmx_inds(1);
        y_displacements = pos_utmy_inds - pos_utmy_inds(1);
    
        % and use that displacement to move the neighborhood and get all
        % indices
        neighborhood_x_coords = repmat(x_inds_init',numel(x_displacements),1) + x_displacements;
        neighborhood_y_coords = repmat(y_inds_init',numel(y_displacements),1) + y_displacements;
    
        % check which neighborhoods are complete in the image
        valid_x = all(neighborhood_x_coords <= size(grid_x,2) & neighborhood_x_coords > 0, 2);
        valid_y = all(neighborhood_y_coords <= size(grid_x,1) & neighborhood_y_coords > 0, 2);
        has_imagery = all([valid_x,valid_y],2);
    
        % take out neighborhoods without imagery
        neighborhood_x_coords = neighborhood_x_coords(has_imagery, :);
        neighborhood_y_coords = neighborhood_y_coords(has_imagery, :);
    
        % just use linear indices for speed and simplicity. As a column vector
        % for now. We can smartly reshape later
        linear_inds = sub2ind(size(grid_x), neighborhood_y_coords(:), neighborhood_x_coords(:));
    
        % make the image into a 2d matrix
        col_image = reshape(satellite_image, [], size(satellite_image,3));
    
        % index the image at the proper location
        col_pix_vals = col_image(linear_inds,:);
    
        % zeros array to fill in pixel values where we have imagery
        pixel_vals = zeros(numel(pos_utmy_inds), size(satellite_image,3));
       
        % now we have to reshape and take the mean to report. Can do all in one line 
        pixel_vals(has_imagery,:) = squeeze(mean(reshape(col_pix_vals,sum(has_imagery), size(neighborhood_x_coords,2), size(satellite_image,3)),2));
    end


    % easiest grabbing of neighborhoods and assessing mean pixel vals is a
    % for loop, but that's slow. Comment this out and try something faster
%     pixel_vals = zeros(numel(pos_utmy_inds), size(satellite_image,3));
%     for i = 1:numel(pos_utmy_inds)
%         % make a circular mask
%         mask = (grid_x - grid_x(1,pos_utmx_inds(i))).^2 + (grid_y - grid_y(pos_utmy_inds(i),1)).^2 <= (neighborhood_size*pix_scale).^2;
%         % grab the pixels and reshape into an array with the number of
%         % columns matching the number of channels
%         mask3d = repmat(mask,1,1,size(satellite_image,3));
%         selected_pix = satellite_image(mask3d);
%         these_pixels = reshape(selected_pix,[],size(satellite_image,3));
% 
%         % collapse via the mean and add to the output array
%         pixel_vals(i,:) = mean(these_pixels);
%     end


end
% function to get the group means of only surface or land pixels

function [mean_utmx,mean_utmy] = loc_mean(utm_x,utm_y,phot_inds)
    % get rid of bottom or noise returns
    good_phots = ismember(phot_inds, [1,3]);
    
    % and calculate the means
    mean_utmx = mean(utm_x(good_phots));
    mean_utmy = mean(utm_y(good_phots));
end