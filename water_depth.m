function [mean_depths, depths_sigma] = water_depth(along_bins,height,photon_ids)
% This function takes refraction-corrected photon heights and uses monte
% carlo simulations to calculate a distribution of water depths within
% designated bins along a track.
%
% IN: 
%
% along_bins: discritization indices for all photons, putting them into
% bins based upon distance along the track. From something like:
% along_bins = discretize(gdats.along{track_ind},bin_edges);
% 
% height: double vector containing the photon return heights for a
% bin on a track of ICESAT-2 data. These should be refraction-corrected
%
% photon_ids: column vector of the same length as input_heights with
% indices designating all classes. Output from cluster_photons.m
%
% OUT: 
%
% mean_depths: mean of depths from distribution of possible water depths
% for each bin.
%
% depths_sigma: standard deviation of depths from distribution of possible
% water depths for each bin.
%
% Written by R. A. Manzuk
% Thursday, January 12, 2023 at 2:06:13 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % start with empty array to get the water depths and sigmas
    mean_depths = zeros(max(along_bins),1);
    depths_sigma = zeros(max(along_bins),1);

    % loop through and get 'em
    for i = 1:numel(mean_depths)
        % just grab these heights
        these_heights = height(along_bins == i);
        these_ids = photon_ids(along_bins == i);

        % get rid of any noise or land pixels
        good_points = these_ids ~= 0 & these_ids ~= 3;

        % only proceed if we have enough other points
        if sum(these_ids == 1) > 1 && sum(these_ids == 2) > 1
            % we want to get the distribution of water depths through monte carlo,
            % so we'll randomly sample a bunch the surface and deep points and subtract
            % them
            % get rid of noisy data
            not_noise = these_ids > 0;
            final_data = these_heights(not_noise);
            final_ids = these_ids(not_noise);

            % define number of samples for monte carlo
            n_samples = 1000;

            all_depths = randsample(final_data(final_ids == 1),n_samples,true) - randsample(final_data(final_ids == 2),n_samples,true);
            % take the mean and standard deviation to report it
            mean_depths(i) = mean(all_depths);
            depths_sigma(i) = std(all_depths);
        else
            mean_depths(i) = 0;
            depths_sigma(i) = 0;
        end
    end
end