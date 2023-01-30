function [centroid_dists, class_centroids] = spectral_clustering(satellite_image, n_clusters)
%
% IN: 
%
% satellite_image: 3d matrix with the multichannel satellite image from
% which you will be calibrating Beer's Law.
%
% n_clusters: number of clusters to be used for the kmeans clustering of
% the satellite imgage
%
% OUT:
%
% centroid_dists: 3D matrix with same number of rows and columns as the
% input satellite image. The third dimension is set by the number of
% clusters. Each 'channel' in this output image is the distance between
% that pixel's spectral data and that particular centroid.
%
% class_centroids: n_clusters x n_channels matrix with the locations of
% each centroid.
%
% Written by R. A. Manzuk
% Monday, January 30, 2023 at 11:31:01 AM
%
% Edited: Monday, January 30, 2023 at 11:31:06 AM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
    % make the satellite into columns
    column_im = reshape(satellite_image, [], size(satellite_image,3));

    % we may have masked out pixels that are 0, so figure out where those
    % are to not include them in clustering
    not_masked = all(column_im > 0, 2);

    % start with empty distance matrix, where we'll leave masked pixels as
    % zero and fill in around them
    distance_cols = zeros(size(column_im,1),n_clusters);


    % cluster the image with k means given the input number of clusters
    [~, class_centroids, ~, distance_cols(not_masked,:)] = kmeans(column_im(not_masked,:), n_clusters);

    % and reshape the distance_columns into the original image shape for
    % export
    centroid_dists = reshape(distance_cols, size(satellite_image,1), size(satellite_image,2), n_clusters);
   
end
