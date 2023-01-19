% just trying to get the hang of things playing with ICESat-2 data and
% relating to Emily's bahamas bathymetry. So some basic blocks to load the
% data, etc.
%
% Created 01/04/2023
% Last edited Friday, January 13, 2023 at 9:34:10 PM
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

%% displaying ICESat tracks on imagery and bathymetry
figure()
subplot(1,2,1)
mapshow(rgb_im,refmat)
hold on
scatter(gdats.utmx{1,1},gdats.utmy{1,1})
scatter(gdats.utmx{1,2},gdats.utmy{1,2})
scatter(gdats.utmx{1,3},gdats.utmy{1,3})
scatter(gdats.utmx{1,4},gdats.utmy{1,4})
scatter(gdats.utmx{1,5},gdats.utmy{1,5})
scatter(gdats.utmx{1,6},gdats.utmy{1,6})
subplot(1,2,2)
% to display on bathymetry, note the tagging of the option to display as a
% surface. mapshow defaults to image inputs.
mapshow(bathy,refmat,'DisplayType','surface')
hold on
scatter(gdats.utmx{1,1},gdats.utmy{1,1})
scatter(gdats.utmx{1,2},gdats.utmy{1,2})
scatter(gdats.utmx{1,3},gdats.utmy{1,3})
scatter(gdats.utmx{1,4},gdats.utmy{1,4})
scatter(gdats.utmx{1,5},gdats.utmy{1,5})
scatter(gdats.utmx{1,6},gdats.utmy{1,6})
a=colorbar;
ylabel(a,'Water depth (m)','FontSize',12,'Rotation',270);
axis tight

%% Extract the bathymetry along one of the track profiles
% first pick a track
track_ind = 1;
% get the starting and ending points in utm space from the track
utmx_start = gdats.utmx{track_ind}(1);
utmx_end = gdats.utmx{track_ind}(end);
utmy_start = gdats.utmy{track_ind}(1);
utmy_end = gdats.utmy{track_ind}(end);

% and use those starting locs to find the closest pixels in the utm grids
[~,utmx_start_closest] = min(abs(utm_x - utmx_start));
[~,utmx_end_closest] = min(abs(utm_x - utmx_end));
[~,utmy_start_closest] = min(abs(utm_y - utmy_start));
[~,utmy_end_closest] = min(abs(utm_y - utmy_end));

% now take the profiles
bathymetry_profile = improfile(bathy,[utmx_start_closest,utmx_end_closest],[utmy_start_closest,utmy_end_closest]);
% and the associated locations with the profile
utmx_profile = improfile(UTM_x,[utmx_start_closest,utmx_end_closest],[utmy_start_closest,utmy_end_closest]);
utmy_profile = improfile(UTM_y,[utmx_start_closest,utmx_end_closest],[utmy_start_closest,utmy_end_closest]);

% get the distance along the profile
d_along_profile = sqrt((utmx_profile - utmx_profile(1)).^2 + (utmy_profile - utmy_profile(1)).^2);
%% distance along the track vs. photon height is where we can do the most work
% first, we'll want to bin things based upon where they fall in distance
% along
% pick a bin size and bin 

bin_size = 200;
bin_edges = [min(gdats.along{track_ind}):bin_size:max(gdats.along{track_ind})];
along_bins = discretize(gdats.along{track_ind},bin_edges);

% use written function decide which photons came from water surface,
% bottom, noise, and mask land
photon_ids = zeros(size(gdats.height{track_ind}));
for i = 1:max(along_bins)
    these_photons = gdats.height{track_ind}(along_bins == i);
    this_mask = gdats.track_mask{track_ind}(along_bins == i);
    photon_ids(along_bins == i) = cluster_photons(these_photons, this_mask, 0.1, 10);
end

% for those that have been identified as being from the bottom, we need to
% correct for refraction
[corrected_height, corrected_utmx, corrected_utmy] = refraction_correct(gdats.height{track_ind}, gdats.utmx{track_ind}, gdats.utmy{track_ind}, gdats.refelev{track_ind}, gdats.azi{track_ind}, photon_ids, along_bins);

% now that we've corrected positions, we need to have a new along track
[adjusted_along] = recalc_along(gdats.along{track_ind}, corrected_utmx, corrected_utmy);

% and redo the edges and centers
bin_edges2 = [min(adjusted_along):bin_size:max(adjusted_along)];
along_bins2 = discretize(adjusted_along,bin_edges2);
bin_centers = (bin_edges2(2:end) + bin_edges2(1:end-1))/2;

% now we just calculate water depth
[mean_depths, depths_sigma] = water_depth(along_bins2,corrected_height,photon_ids); 

