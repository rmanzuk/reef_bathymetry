% As I try to contain all elements needed for bathymetry analysis of a
% shallow platform via ICESat-2 into a single class, I'm using the Bahamas
% as a test case. This brief script runs through basic class construction
% and calculation functions
%
% Ryan A. Manzuk 09/08/2023
%%
% make up a platoform for the bahamas data
bahamas = icesat_platform;
bahamas.nameid = 'Bahamas';
bahamas.track_dir = '/Users/ryan/Dropbox (Princeton)/modern_reef_bathymetry/platforms/bahamas/icesat_tracks';
bahamas.sat_im_path = '/Users/ryan/Dropbox (Princeton)/modern_reef_bathymetry/platforms/bahamas/sat_im.tif';
bahamas.sat_mask_path = '/Users/ryan/Dropbox (Princeton)/modern_reef_bathymetry/platforms/bahamas/mask.tif';

% use functions to load and adjust properties
bahamas.load_sat_im;
bahamas.load_sat_mask;
bahamas.make_sat_rgb;
bahamas.load_tracks;
bahamas.set_track_utm;

% cluster the satellite image
bahamas.cluster_sat;

% apply satellite image mask to tracks
bahamas.set_track_mask;

% calculate water depth along all tracks
bahamas.calc_depth;

% get the pixel values and cluster assignments from the track locations
% where we took water depths
bahamas.pixel_vals('rawsat')
bahamas.pixel_vals('cluster');

% TO DO: 
% 1) make the calc_depth function more readable and modular
% 2) Download more ICESat data
% 3) Deriving tidal structure from ICESat + sat photos?