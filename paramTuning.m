% After creating a baseline method for water depth from ICESAT-2 photon
% returns, need to tune some of the parameters for best results
%
% Created 01/12/2023
% Last edited Tuesday, January 17, 2023 at 11:01:01 AM
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

%% main crux is properly clustering the photons
% first pick a track
track_ind = 1;

% from manual inspection of this track, I know a really important section
% of this track for clustering occurs between approx 2.04e5 and 2.0e5
% meters along. this section contains a transition from land to shallow
% water, which is the hardest to accurately cluster
within_sect = gdats.along{track_ind} > 2.04e5 & gdats.along{track_ind} < 2.07e5;
these_along = gdats.along{track_ind}(within_sect);
these_heights = gdats.height{track_ind}(within_sect);
this_mask = gdats.track_mask{track_ind}(within_sect);

% four tunable parameters, but let's just fix the neighbor distance for
% now and manipulate bin sizes, and n_neighbors
dist_thresh = 0.1;
bin_sizes = [10,15,20,50,100,200,500];
frac_neighbors = [0.01,0.025,0.05,0.1,0.15,0.2];

% loop through all the parameter spaces and get some results
photon_ids = zeros(numel(these_along), numel(bin_sizes), numel(n_neighbors));
for i = 1:numel(bin_sizes)
    bin_size = bin_sizes(i);
    bin_edges = [min(these_along):bin_size:max(these_along)];
    along_bins = discretize(these_along,bin_edges);
    for j = 1:numel(frac_neighbors)
        frac_thresh = frac_neighbors(j);
        % use written function decide which photons came from water surface,
        % bottom, land, etc.
        for k = 1:max(along_bins)
            these_photons = these_heights(along_bins == k);
            bin_mask = this_mask(along_bins == k);
            photon_ids(along_bins == k,i,j) = cluster_photons(these_photons, bin_mask, dist_thresh, frac_thresh);
        end
    end
end
%% plot up what those clusterings look like
figure()
for i = 1:numel(bin_sizes)
    for j = 1:numel(frac_neighbors)
        subplot(numel(bin_sizes), numel(n_neighbors), (numel(n_neighbors)*(i-1))+j)
        hold on
        scatter(these_along(photon_ids(:,i,j) == 1),these_heights(photon_ids(:,i,j) == 1),5,'filled')
        scatter(these_along(photon_ids(:,i,j) == 2),these_heights(photon_ids(:,i,j) == 2),5,'filled')
        scatter(these_along(photon_ids(:,i,j) == 3),these_heights(photon_ids(:,i,j) == 3),5,'filled')
        scatter(these_along(photon_ids(:,i,j) == 0),these_heights(photon_ids(:,i,j) == 0),5,'filled')
        hold off
        ylim([-34 -25])
        title([num2str(bin_sizes(i)) 'm bins, ' num2str(frac_neighbors(j)) ' pts'], 'FontSize', 10);
    end
end