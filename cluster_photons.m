function [photon_ids] = cluster_photons(input_heights, input_mask, dist_thresh, frac_points)
% This function looks within an along-track length bin of ICESAT-2 data and
% distinguishes which photon returns are noise, land, sea surface, and sea
% bottom.
%
% IN: 
%
% input_heights: double vector containing the photon return heights for a
% bin on a track of ICESAT-2 data. These should not be corrected in any
% way.
%
% input_mask: logical vector of same length as input_heights containing the
% corresponding land mask indication for all photon returns.
% 
% dist_thresh: Distance within which to search for neighbors when
% determining which points are noise (in meters). 
%
% frac_points: Threshold for minimum fraction of points in the bin that
% must be within the dist_thresh to be consideren neighbors.
%
% OUT: 
%
% photon_ids: column vector of the same length as input_heights with
% indices designating all classes: 
% 0: noise
% 1: sea surface
% 2: sea bottom
% 3: land
%
% Written by R. A. Manzuk
% Thursday, January 12, 2023 at 11:44:41 AM
%
% Edited: Tuesday, January 17, 2023 at 11:01:22 AM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % make an empty array for the inds
    photon_ids = zeros(size(input_heights));

    % put in the land indices right away
    photon_ids(~input_mask) = 3;

    % and store these overall sea inds for later
    sea_inds = find(input_mask == true);

    % and lets remove land points from further consideration
    sea_heights = input_heights(input_mask);

    % pairwise distances between points is a good place to start
    dists = pdist2(sea_heights, sea_heights);

    % how many distances are less than the threshold for each point
    num_neighbors = sum(dists < dist_thresh, 2);

    % make a logical to id noise and not_noise
    n_points = round(frac_points * numel(input_heights));
    not_noise = num_neighbors > n_points;  

    % get rid of noisy data
    denoised = sea_heights(not_noise);
    
    % we'll only proceed to clustering bottom and surface if there
    % are any non-noise points left
    if numel(denoised) > 2
        % for now let's id surface and bottom via kmeans
        [clusters, centroids] = kmeans(denoised,2);

        % the surface cluster is that which has the higher centroid
        [~,surf_ind] = max(centroids);
        
        % we always want surf_ind to be 1, so if it's not, we need to
        % switch the cluster assignments
        if surf_ind ~= 1
            % grab the opposites as surface
            surfies = clusters == 2;
            % reassign
            clusters(surfies) = 1;
            clusters(~surfies) = 2;
        end

        % now we put them into photon_ids
        photon_ids(sea_inds(not_noise)) = clusters;
    end
end
% code for distinguishing land used prior to implementation of land masks
% on Friday, January 13, 2023 at 9:05:52 PM
%     %dig into the distances to centroids to decide if this was land
%     % start by taking the ratio of the distance to the further centroid/the
%     % distance to the nearer centroid
%     dist_ratios = [cent_dists(clusters == 1,2)./cent_dists(clusters == 1,1);...
%         cent_dists(clusters == 2,1)./cent_dists(clusters == 2,2)];
%     % land should have some photons where the distance ratio approaches
%     % one, as in they are equidistant from each centroid
%     % right now we'll just set a threshold based upone the minimum of the
%     % distance ratios, maybe dig more into the distribution of distance
%     % ratios in the future
%     min_thresh = 1.1;
%     % if we are less than this thresh, we can set all the inds of non-noise
%     % to land (3)
%     if min(dist_ratios) < min_thresh
%         photon_ids(not_noise) = 3;
%     else