% Now that I've gotten used to dealing with tracks and have a minimally
% functional method for water depth, time to run through thngs and
% calibrate Beer's Law.
%
% Created 01/17/2023
% Last edited Tuesday, January 17, 2023 at 11:00:49 AM
% R. A. Manzuk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading of bathymetry data from StartUpScript included in the EmilyDataCube

% How to load data and prep matrices.
load('/Users/ryan/Dropbox (Princeton)/modern_reef_bathymetry/IceSat2/Bahamas/EmilyDataCube/data_cube_v2.mat') % variables: 'bathy','mask','utm_x','utm_y'

% If you want the UTM coordinates in a 2D grid:
[UTM_x,UTM_y] = meshgrid(utm_x,utm_y);

% optional: if you want the imagery (already sampled to the same grid as
% bathymetry, mask, and UTM coordinates). Note that this is big (~1.2 Gb)
[RE, refmat, bbox] = geotiffread('/Users/ryan/Dropbox (Princeton)/modern_reef_bathymetry/IceSat2/Bahamas/EmilyDataCube/RE_Clipped.tif');

% good to just make a quick rgb version of the satelite image for display
rgb_im = imlocalbrighten(RE(:,:,[3,2,1]));
%% loading an ICESat track
%
load('/Users/ryan/Dropbox (Princeton)/modern_reef_bathymetry/IceSat2/Bahamas/gdats-ATL03_20181103075345_05450107_003_01.mat');

% first thing to note is utmx, utmy fields of gdats appear to be missing,
% so fill them in based upon a conversion from the lat long field
for i = 1:numel(gdats.lon)
    % get the utm zone from the coordinates
    positions = [gdats.lat{i},gdats.lon{i}];
    utm_zone = utmzone(positions);
    % get the geoid of the zon and construct its projection structure
    [ellipsoid,estr] = utmgeoid(utm_zone);
    utmstruct = defaultm('utm');
    utmstruct.zone = utm_zone;
    utmstruct.geoid = ellipsoid;
    utmstruct = defaultm(utmstruct);
    % and just do the conversion
    [gdats.utmx{i},gdats.utmy{i}] = mfwdtran(utmstruct,gdats.lat{i},gdats.lon{i});

    % while we're at it, let's mask our land returns
    gdats.track_mask{i} = track_land_mask(gdats.utmx{i}, gdats.utmy{i}, mask, UTM_x, UTM_y);
end

%% first step is getting all the water depths for each track
% first, we'll want to bin things based upon where they fall in distance
% along
% pick a bin size and bin 
bin_size = 50;

% Set the algorithm thresholds learned from tuning
dist_thresh = 0.1;
frac_neighbors = 0.025;

% loop through and do this for all tracks
for i = 1:numel(gdats.lon)
    % i is the track_ind
    track_ind = i;
    % bin that thing
    bin_edges = [min(gdats.along{track_ind}):bin_size:max(gdats.along{track_ind})];
    along_bins = discretize(gdats.along{track_ind},bin_edges);
    % use written function decide which photons came from water surface,
    % bottom, noise, and mask land
    photon_ids = zeros(size(gdats.height{track_ind}));
    for j = 1:max(along_bins)
        these_photons = gdats.height{track_ind}(along_bins == j);
        this_mask = gdats.track_mask{track_ind}(along_bins == j);
        photon_ids(along_bins == j) = cluster_photons(these_photons, this_mask, dist_thresh, frac_neighbors);
    end
    
    % place the photon_ids into the struct
    gdats.photon_ids{track_ind} = photon_ids;

    % for those that have been identified as being from the bottom, we need to
    % correct for refraction
    [corrected_height, corrected_utmx, corrected_utmy] = refraction_correct(gdats.height{track_ind}, gdats.utmx{track_ind}, gdats.utmy{track_ind}, gdats.refelev{track_ind}, gdats.azi{track_ind}, photon_ids, along_bins);
    
    % the corrected things into the struct
    gdats.corrected_height{track_ind} = corrected_height;
    gdats.corrected_utmx{track_ind} = corrected_utmx;
    gdats.corrected_utmy{track_ind} = corrected_utmy;

    % now that we've corrected positions, we need to have a new along track
    [adjusted_along] = recalc_along(gdats.along{track_ind}, corrected_utmx, corrected_utmy);

    % put that in the struct too
    gdats.adjusted_along{track_ind} = adjusted_along;
    
    % and redo the edges and centers
    bin_edges2 = [min(adjusted_along):bin_size:max(adjusted_along)];
    along_bins2 = discretize(adjusted_along,bin_edges2);
    bin_centers = (bin_edges2(2:end) + bin_edges2(1:end-1))/2;
    
    % now we just calculate water depth
    [mean_depth, depth_sigma] = water_depth(along_bins2,corrected_height,photon_ids);

    % place the depths, std devs and the associated bin centers into the struct 
    gdats.mean_depth{track_ind} = mean_depth;
    gdats.depth_sigma{track_ind} = depth_sigma;
    gdats.bin_centers{track_ind} = bin_centers;
end

%% now that we have water depths for all tracks need to get their associated pixel values from the satellite image
% select a neighboorhood size for pixel averaging in the image
neighborhood_size = 3;

% loop through and do this for all tracks
for i = 1:numel(gdats.lon)
        % i is the track_ind
    track_ind = i;
    % bin that thing
    bin_edges = [min(gdats.adjusted_along{track_ind}):bin_size:max(gdats.adjusted_along{track_ind})];
    along_bins = discretize(gdats.adjusted_along{track_ind},bin_edges);

    % extract pixel values
    [gdats.pixel_vals{track_ind}] = track_pixels(along_bins, gdats.photon_ids{track_ind},  gdats.utmx{track_ind},  gdats.utmy{track_ind}, RE, UTM_x, UTM_y, neighborhood_size);
end

% get all tracks in one place for good measure
all_depths = [gdats.mean_depth{1}; gdats.mean_depth{2}; gdats.mean_depth{3};...
    gdats.mean_depth{4}; gdats.mean_depth{5}; gdats.mean_depth{6}];
all_pixels = [gdats.pixel_vals{1}; gdats.pixel_vals{2}; gdats.pixel_vals{3};...
    gdats.pixel_vals{4}; gdats.pixel_vals{5}; gdats.pixel_vals{6}];
cleaned_depths = all_depths(all_depths > 0 & all_depths < 10 & all_pixels(:,1) > 0);
cleaned_pixels = all_pixels(all_depths > 0 & all_depths < 10 & all_pixels(:,1) > 0, :);

%% a nice plot of water depth vs color values


figure()
subplot(1,2,1)
scatter(cleaned_pixels(:,3), cleaned_depths, 10, 'filled', 'DisplayName', 'Band 3');
hold on
scatter(cleaned_pixels(:,4), cleaned_depths, 10, 'filled', 'DisplayName', 'Band 4');
scatter(cleaned_pixels(:,5), cleaned_depths, 10, 'filled', 'DisplayName', 'Band 5');
hold off
xlabel('Pixel value')
ylabel('Water Depth [m]')
ylim([0 10])
legend()
%  and include the log version
subplot(1,2,2)
scatter(log10(cleaned_pixels(:,3)), cleaned_depths, 10, 'filled', 'DisplayName', 'Band 3');
hold on
scatter(log10(cleaned_pixels(:,4)), cleaned_depths, 10, 'filled', 'DisplayName', 'Band 4');
scatter(log10(cleaned_pixels(:,5)), cleaned_depths, 10, 'filled', 'DisplayName', 'Band 5');
hold off
xlabel('log(Pixel value)')
ylabel('Water Depth [m]')
ylim([0 10])

%% band ratios are better for correcting for bottom types when assessing depth from pixel value

n_channels = size(RE,3);

% we'll take all ratios and perform a sensitivity analysis
combos = nchoosek(1:n_channels,2);
permutations = reshape(combos(:,perms(1:2)),[],2);

% note its the log band ratios we take
band_ratios = log10(cleaned_pixels(:,permutations(:,1)))./log10(cleaned_pixels(:,permutations(:,2)));

[r_squareds, best_ratios] = variance_reduction(band_ratios,cleaned_depths);

% and make a plot that helps us visualize the combination of bands that add
% the most to the regression
figure()
hold on
scatter(1:size(band_ratios,2),max(r_squareds), 30, 'filled')
text([1:size(band_ratios,2)] - 0.25, max(r_squareds) + 0.01, [num2str(permutations(best_ratios,1)) repmat('/',size(band_ratios,2),1) num2str(permutations(best_ratios,2))])
ylabel('R^2')

