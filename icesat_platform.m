classdef icesat_platform < handle
   properties
       % what to call this platform
       nameid 
       % mostly a lot of paths to files as properties, like where track
       % data are stored
       track_dir
       % the satellite image path
       sat_im_path
       % a mask for the satellite image
       sat_mask_path
       % a place to store the satellite image and mask once loaded in
       sat_im
       sat_mask
       % grids for the utm x and y locations of all pixel centroids
       utmx_grid
       utmy_grid
       % the scale of pixels in the satellite image
       sat_pixel_scale
       % a place to store icesat tracks once loaded in
       gdats
       % just an RGB version of the satellite image for display, etc.
       sat_rgb
       % cluster assignments for pixels in the satellite image
       sat_clusters
   end
   %%
   methods
        %%
        % load in the satellite image and get geographic data
        function load_sat_im(obj)

            % read in the geotiff, storing spatial referencing info
            [obj.sat_im, spat_ref] = readgeoraster(obj.sat_im_path);

            % use the spatial referencing to get a utm grid for the pixels
            edges_utmx = spat_ref.XWorldLimits(1):spat_ref.CellExtentInWorldX:spat_ref.XWorldLimits(2);
            edges_utmy = spat_ref.YWorldLimits(1):spat_ref.CellExtentInWorldY:spat_ref.YWorldLimits(2);
            
            % average the pixel edges to get centroids
            centers_utmx = (edges_utmx(2:end) + edges_utmx(1:end-1))./2;
            centers_utmy = (edges_utmy(2:end) + edges_utmy(1:end-1))./2;

            % make into meshgrid for export
            [obj.utmx_grid, obj.utmy_grid] = meshgrid(centers_utmx, flip(centers_utmy));
            
            % and export the scale
            obj.sat_pixel_scale = (spat_ref.CellExtentInWorldX + spat_ref.CellExtentInWorldY)/2;
        end
        
        %%
        % simple to load in the satellite mask
        function load_sat_mask(obj)
            obj.sat_mask = im2double(imread(obj.sat_mask_path));
        end
        
        %%
        % for display or other purposes, can be good to have an rgb image
        % from the five channel rapideye image
        function make_sat_rgb(obj)
            % if we haven't loaded the sat image, do so
            if isempty(obj.sat_im)
                obj.load_sat_im;
            end

            % and we're good to go. We can contrast shift this one too
            obj.sat_rgb = imlocalbrighten(obj.sat_im(:,:,[3,2,1]));
        end

        %%
        % load in all of the icesat tracks available into one nice struct
        function load_tracks(obj)
            % make a directory for any .mat file in the track folder
            track_directory = dir([obj.track_dir '/*.mat']);
            
            % load in the first one to start the structure and get ready to
            % track fields
            raw_gdats = load(fullfile(obj.track_dir, track_directory(1).name));
            raw_gdats = raw_gdats.gdats;
            fields = fieldnames(raw_gdats);

            % if there are any empty tracks, we want to ignore them
            empty_track = cellfun("isempty", raw_gdats.(fields{1}));

            % and get the object gdats started, clearing empty fields if
            % they are there
            fields_adj = fields(:).';
            fields_adj{2,1} = {};
            obj.gdats = struct(fields_adj{:});
            for i = 1:numel(fields)
                obj.gdats(1).(fields{i}) = raw_gdats.(fields{i})(~empty_track);
            end
            
            % we want to id which pass each track came from, just with a
            % numerical identifier
            obj.gdats.passid = ones(1,sum(~empty_track));
      
            % and loop through everything
            for i = 2:numel(track_directory)
                new_gdats = load(fullfile(obj.track_dir, track_directory(i).name));
                new_gdats = new_gdats.gdats;
                
                % some tracks might be empty, so we need to be able to
                % ignore them
                empty_track = cellfun("isempty", new_gdats.(fields{1}));

                % gotta loop through all the fields too 
                for j = 1:numel(fields)
                    obj.gdats.(fields{j}) = [obj.gdats.(fields{j}), new_gdats.(fields{j})(~empty_track)];
                end

                % put in the passid
                obj.gdats.passid = [obj.gdats.passid, i*ones(1,sum(~empty_track))];
            end
        end

        %%
        % it's possible that the utmx and y fields of the struct are empty,
        % so we can set them
        function set_track_utm(obj)
            % just going to loop through and do this for all
            for i = 1:numel(obj.gdats.lon)
                % get the utm zone from the coordinates
                positions = [obj.gdats.lat{i},obj.gdats.lon{i}];
                utm_zone = utmzone(positions);
                % get the geoid of the zon and construct its projection structure
                [ellipsoid,~] = utmgeoid(utm_zone);
                utmstruct = defaultm('utm');
                utmstruct.zone = utm_zone;
                utmstruct.geoid = ellipsoid;
                utmstruct = defaultm(utmstruct);
                % and just do the conversion
                [obj.gdats.utmx{i},obj.gdats.utmy{i}] = projfwd(utmstruct,obj.gdats.lat{i},obj.gdats.lon{i});
            end
        end
        %%
        % apply the mask to each track so we know wich photons we don't
        % really need to worry about
        function set_track_mask(obj)
            % check that we loaded in the mask
            if isempty(obj.sat_mask)
                obj.load_sat_mask;
            end

            % and loop through each track and use a function to apply the
            % mask
            for i = 1:numel(obj.gdats.lon)
                obj.gdats.track_mask{i} = track_land_mask(obj.gdats.utmx{i}, obj.gdats.utmy{i}, obj.sat_mask, obj.utmx_grid, obj.utmy_grid);
            end
        end
        
        %%
        % cluster the satellite image
        function cluster_sat(obj)
            % hard code number of clusters for now
            n_clusters = 5;

            % make sure we have the satellite image and mask loaded
            if isempty(obj.sat_im)
                obj.load_sat_im
            end

            if isempty(obj.sat_mask)
                obj.load_sat_mask
            end
            
            % don't want land included in the clustering
            sat_masked = im2double(obj.sat_im) .* obj.sat_mask;
            
            % run function to get clustering results
            [centroid_dists, ~] = spectral_clustering(sat_masked, n_clusters);
            
            % and the hard cluster assignments are just the indices of the minima of
            % the centroid distances.
            [~, obj.sat_clusters] = min(centroid_dists,[],3);
        end
        %%
        % get the water depths for each track
        function calc_depth(obj)

            % only proceed if we have tracks loaded in. If we don't load
            % them
            if isempty(obj.gdats)
                obj.load_tracks;
            end
           
            % and we need to have utm locations, so double check for those
            if isempty(obj.gdats.utmx)
                obj.set_track_utm;
            end

            % hard code in a few parameters for right now
            bin_size = 50;
            dist_thresh = 0.1;
            frac_neighbors = 0.025;

            % ready to go....we're gonna loop through all the tracks and
            % just cluster the photons first
            for i = 1:numel(obj.gdats.lon)

                % i is the track_ind
                track_ind = i;

                % bin that thing
                bin_edges = [min(obj.gdats.along{track_ind}):bin_size:max(obj.gdats.along{track_ind})];
                along_bins = discretize(obj.gdats.along{track_ind},bin_edges);

                % use written function decide which photons came from water surface,
                % bottom, noise, and mask land
                photon_ids = zeros(size(obj.gdats.height{track_ind}));
                for j = 1:max(along_bins)
                    these_photons = obj.gdats.height{track_ind}(along_bins == j);
                    this_mask = obj.gdats.track_mask{track_ind}(along_bins == j);
                    photon_ids(along_bins == j) = cluster_photons(these_photons, this_mask, dist_thresh, frac_neighbors);
                end
                
                % place the photon_ids into the struct
                obj.gdats.photon_ids{track_ind} = photon_ids;
            end

            % and now we can 'tidal' correct the heights by fitting a plain
            % to the surface-identified photons and removing its tilt
            obj.tidal_correct;

            % and loop back through the tracks to correct for refraction
            % and calculate depth
            for i = 1:numel(obj.gdats.lon)

                % i is the track_ind
                track_ind = i;

                % bin that thing
                bin_edges = [min(obj.gdats.along{track_ind}):bin_size:max(obj.gdats.along{track_ind})];
                along_bins = discretize(obj.gdats.along{track_ind},bin_edges);
            
                % for those that have been identified as being from the bottom, we need to
                % correct for refraction
                [corrected_height, corrected_utmx, corrected_utmy] = refraction_correct(obj.gdats.desloped_heights{track_ind}, obj.gdats.utmx{track_ind}, obj.gdats.utmy{track_ind}, obj.gdats.refelev{track_ind}, obj.gdats.azi{track_ind}, obj.gdats.photon_ids{track_ind}, along_bins);
                
                % put the corrected things into the struct
                obj.gdats.corrected_height{track_ind} = corrected_height;
                obj.gdats.corrected_utmx{track_ind} = corrected_utmx;
                obj.gdats.corrected_utmy{track_ind} = corrected_utmy;
            
                % now that we've corrected positions, we need to have a new along track
                [adjusted_along] = recalc_along(obj.gdats.along{track_ind}, corrected_utmx, corrected_utmy);
            
                % put that in the struct too
                obj.gdats.adjusted_along{track_ind} = adjusted_along;
                
                % and redo the edges and centers
                bin_edges2 = min(adjusted_along):bin_size:max(adjusted_along);
                along_bins2 = discretize(adjusted_along,bin_edges2);
                bin_centers = (bin_edges2(2:end) + bin_edges2(1:end-1))/2;
                
                % now we just calculate water depth
                [mean_depth, depth_sigma] = water_depth(along_bins2,corrected_height,obj.gdats.photon_ids{track_ind});
            
                % place the depths, std devs and the associated bin centers into the struct 
                obj.gdats.mean_depth{track_ind} = mean_depth;
                obj.gdats.depth_sigma{track_ind} = depth_sigma;
                obj.gdats.bin_centers{track_ind} = bin_centers;
            end
        end
   
        %% 
        % put the pixel values from the satellite image for each associated
        % track location water depth
        function pixel_vals(obj, source)
            % user can input the source as either the raw satellite image
            % or the cluster map
            if strcmp(source, 'rawsat')
                % check for the satellite im
                if isempty(obj.sat_im)
                    obj.load_sat_im
                end
            elseif strcmp(source, 'cluster')
                % check for the cluster im
                if isempty(obj.sat_clusters)
                    obj.cluster_sat
                end
            else
                error('invalid source type')
            end

            % ready to go. Start with hard code of neighborhood size and
            % bin size. bin size needs to be the same as from calc_depth.
            neighborhood_size = 3;
            bin_size = 50;

            % loop through and do this for all tracks
            for i = 1:numel(obj.gdats.lon)
                % i is the track_ind
                track_ind = i;
                % bin that thing
                bin_edges = [min(obj.gdats.adjusted_along{track_ind}):bin_size:max(obj.gdats.adjusted_along{track_ind})];
                along_bins = discretize(obj.gdats.adjusted_along{track_ind},bin_edges);
                
                % if user wanted pixel values, extract pixel values
                if strcmp(source, 'rawsat')
                    [obj.gdats.pixel_vals{track_ind}] = track_pixels(along_bins, obj.gdats.photon_ids{track_ind},  obj.gdats.utmx{track_ind},  obj.gdats.utmy{track_ind}, obj.sat_im, obj.utmx_grid, obj.utmy_grid, neighborhood_size);
                end

                % if user wanted cluster assignments, extract cluster
                % assignments, but neighborhood will be 0
                if strcmp(source, 'cluster')
                    [obj.gdats.centroid_dists{track_ind}] = track_pixels(along_bins, obj.gdats.photon_ids{track_ind},  obj.gdats.utmx{track_ind},  obj.gdats.utmy{track_ind}, obj.sat_clusters, obj.utmx_grid, obj.utmy_grid, 0);
                end
            end
        end
    
        %%
        % display all tracks on the satellite imgage
        function plot_tracks(obj)

            % first make sure we have the rgb satellite image and the
            % tracks
            if isempty(obj.gdats)
                obj.load_tracks;
            end
            
            if isempty(obj.sat_rgb)
                obj.make_sat_rgb;
            end

            % the upper left corner of the image will be the origin for
            % plotting tracks, so extract that from the grids
            origin_x = obj.utmx_grid(1,1);
            origin_y = obj.utmy_grid(1,1);

            % and we're ready to go
            figure()

            % show the image
            imshow(obj.sat_rgb)
            hold on
            
            % loop through all the tracks
            for i = 1:numel(obj.gdats.lon)
                
                % extract the coordinates and subtract the origin
                this_x = obj.gdats.utmx{i} - origin_x;
                this_y = -(obj.gdats.utmy{i} - origin_y);
                
                % in stead of scattering all the points, we can just plot a
                % line of best fit
                ones_array = ones(numel(this_x),1);
                params = [this_x, ones_array]\this_y;

                plot([this_x(1), this_x(end)]./obj.sat_pixel_scale, [this_x(1)*params(1) + params(2), this_x(end)*params(1) + params(2)]./obj.sat_pixel_scale, 'Color',[1, 1, 1])
            end

        end
        
        %%
        % fit a plane to the surface photons for each set of tracks and
        % snap them back to having a horizontal surface
        function tidal_correct(obj)

            % make sure we have track data, otherwise load it
            if isempty(obj.gdats)
                obj.load_tracks
            end

            % we will loop through all passids, fitting once per passid
            for i = 1:max(unique(obj.gdats.passid))
                
                % get x,y,z data for the photons on this pass, and their
                % surface vs. bottom vs. land id
                pass_heights = vertcat(obj.gdats.height{obj.gdats.passid == i});
                pass_utmx = vertcat(obj.gdats.utmx{obj.gdats.passid == i});
                pass_utmy = vertcat(obj.gdats.utmy{obj.gdats.passid == i});
                pass_ids = vertcat(obj.gdats.photon_ids{obj.gdats.passid == i});

                % get just the surface ones
                surface_heights = pass_heights(pass_ids == 1);
                surface_utmx = pass_utmx(pass_ids == 1);
                surface_utmy = pass_utmy(pass_ids == 1);

                % fit a plane with least squares
                ones_vec = ones(numel(surface_utmy), 1);
                params = [surface_utmx, surface_utmy, ones_vec]\surface_heights;

                % we'll just use the mean height of the plane as where
                % we'll pull them all back to
                mean_height = mean(surface_utmx*params(1) + surface_utmy*params(2) + params(3));
                correction_factors = mean_height - (surface_utmx*params(1) + surface_utmy*params(2) + params(3));

                % now we just need to put them back in the struct
                desloped_pass_heights = pass_heights;
                desloped_pass_heights(pass_ids == 1) = pass_heights(pass_ids == 1) + correction_factors;
                
                % need to do some gymnastics to put them back into a cell
                % array of the same size as before
                array_lengths = cellfun(@numel, obj.gdats.height(obj.gdats.passid == i), 'UniformOutput',true);
                desloped_heights_cell = {};
                for j = 1:numel(array_lengths)
                    if j == 1
                        desloped_heights_cell{j} = desloped_pass_heights(1:array_lengths(j));
                    else
                        desloped_heights_cell{j} = desloped_pass_heights((sum(array_lengths(1:j-1))+1):sum(array_lengths(1:j)));
                    end
                end

                obj.gdats.desloped_heights(obj.gdats.passid == i) = desloped_heights_cell;

            end
        end
   end
end