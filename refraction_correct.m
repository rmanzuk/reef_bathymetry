function [corrected_height, corrected_utmx, corrected_utmy] = refraction_correct(height, utmx, utmy, refelev, azi, photon_ids, along_bins)
% This function takes photons (already with land, surface, etc.
% distinctions) and geometrically corrects for distortions in sea bottom
% photon returns that stem from the higher refractive index of water. Based
% largely upon: 
% Parrish, Christopher E., et al. "Validation of ICESat-2 ATLAS bathymetry
% and analysis of ATLAS’s bathymetric mapping performance." Remote sensing
% 11.14 (2019): 1634.
%
% IN: 
%
% height: double vector containing the photon return heights for a
% bin on a track of ICESAT-2 data. These should not be corrected in any
% way.
% 
% utmx: double vector with the easting coordinates in UTM for photon
% returns.
%
% utmy: double vector with the northing coordinates in UTM for photon
% returns.
%
% refelev: double vector with the reference elevation (heading) angle for
% photon returns in radians.
%
% azi: double vector with the azimuth for the laser (relative to UTM grid)
% in radians
%
% photon_ids: column vector of the same length as input_heights with
% indices designating all classes. Output from cluster_photons.m
%
% along_bins: discritization indices for all photons, putting them into
% bins based upon distance along the track. From something like:
% along_bins = discretize(gdats.along{track_ind},bin_edges);
%
% OUT: 
%
% corrected_height: double vector with heights for all photon returns,
% where all bottom grouped photons have been corrected for refraction.
%
% corrected_utmx: double vector with easting coordinates for all photon
% returns, where all bottom grouped photons have been corrected for
% refraction.
% 
% corrected_utmy: double vector with northing coordinates for all photon
% returns, where all bottom grouped photons have been corrected for
% refraction.
% 
%
% Written by R. A. Manzuk
% Thursday, January 12, 2023 at 12:02:06 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    % we can start by defining our indices of refraction
    n_sea = 1.34116;
    n_air = 1.00029;

    % we need the angle of incidence, which Parrish defines as: 
    % pi/2 - ref_elev
    theta1 = pi/2 - refelev(photon_ids == 2);

    % from the angle of incidence, we can use Snell's Law to calculate the
    % angle at which the photons are traveling under water 
    theta2 = asin((n_air*sin(theta1))/n_sea);

    % we also need a recorded depth for each photon, so we'll do that
    % roughly by going within bins and getting a surface height
    % make an empty vector for the depth of all bottom photons
    initial_depth = [];
    for i = 1:max(along_bins)
        these_photons = photon_ids(along_bins == i);

        % first just ask if there are any bottom photons here to worry
        % about
        if sum(these_photons == 2) > 0 
            % if so, we need to know the surface elevation within this bin,
            % we can just take the mean for now, think about something more
            % nuanced later
            bin_heights = height(along_bins == i);
            mean_surf = mean(bin_heights(these_photons == 1));

            % so now we can get the initial depths
            initial_depth = [initial_depth; (mean_surf - bin_heights(these_photons == 2))];
        end
    end

    % time to proceed with the correction
    % find the uncorrected distance traveled by the photon underwater (hypotenuse)
    S = initial_depth ./ cos(theta1);

    % Snell's law to find the actual distance traveled
    R = S * (n_air/n_sea);
    
    % need the angle ratio
    phi = theta1 - theta2;

    % for geometry at the seafloor, need the complementary angle to theta1
    lambda = (pi/2) - theta1;

    % big step, law of cosines to find the distance between the apparent
    % position of the photon and its actual position.
    P = sqrt(R.^2 + S.^2 - 2*R.*S.*cos(theta1-theta2));

    % missing two more angles, because beta will be essential for the
    % correction.
    alpha = asin((R.*sin(phi))./P);
    beta = lambda - alpha;

    % finally, we can correct!
    % depth correction is easy
    delta_Z = P.*sin(beta);

    % to get the easting and northing corrections, first start with the
    % along-track correction
    delta_y = P.*cos(beta);
    % which  can be related to easting and northing from the azimuth of the
    % laser
    delta_E = delta_y.*sin(azi(photon_ids ==2));
    delta_N = delta_y.*cos(azi(photon_ids ==2));

    % just apply the corrections
    corrected_height = height;
    corrected_height(photon_ids == 2) = height(photon_ids == 2) + delta_Z;

    corrected_utmx = utmx;
    corrected_utmx(photon_ids == 2) = utmx(photon_ids == 2) + delta_E;

    corrected_utmy = utmy;
    corrected_utmy(photon_ids == 2) = utmy(photon_ids == 2) + delta_N;
end