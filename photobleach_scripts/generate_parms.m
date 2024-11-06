function generate_parms(cy3_excite, cy5_excite, exposure_time, delay_time, starting_frame_number,end_of_movie, output_folder)
    % Function to generate parameter settings and save them to a file.
    % Inputs:
    %   - cy3_excite: Number of Cy3 excitation frames
    %   - cy5_excite: Number of Cy5 excitation frames
    %   - exposure_time: Exposure time for each frame
    %   - delay_time: Delay time between frames
    %   - starting_frame_number: frame number defining start of photobleach
    %   analysis
    %   - output_folder: Folder to save the output parameters
    
    % Check if the output folder exists; if not, create it
    if ~isfolder(output_folder)
        mkdir(output_folder)
    end
    
    % Set the output file path
    dataout = fullfile(pwd, output_folder);
    
    % Define parameters for image acquisition and analysis
    number_cy3_excitation_frames = cy3_excite;       % Number of Cy3 excitation frames
    number_cy5_excitation_frames = cy5_excite;       % Number of Cy5 excitation frames
    exposuret = exposure_time;                       % Exposure time per frame
    delay = delay_time;                              % Delay time between exposures
    starting_frame = starting_frame_number;          % Frame number to start analyzing movie  
    last_movie_frame = end_of_movie;                 % Frame number to end analyzing movie 
    border = 20;                                     % Width of the boundary between channels (in pixels), default = 20
    boxwidth = 4;                                    % Width of the box around points (in pixels), default = 4
    rectangle = 8;                                   % Box size for reference point least squares fit (LSF) transformation, default = 8
    HighCutoff = 99.5;                                 % Upper limit for spot intensity and size distribution, default = 99
    LowCutoff = 0.05;                               % Lower limit for spot intensity and size distribution, default = 0.05
    drkwindow = 10;                                  % Window size for observing nonspecific binding events (in pixels), default = 10
    area = 12;                                       % Area of the box defining the spot (12x12 pixels), default = 12
    minimum_distance_between_centers = 6;            % Minimum distance between spot centers (in pixels), default = 6
    BGcut = 80;                                      % Percentage of average field of view (FOV) intensity for identifying local maxima, default = 80
    Gauss_Drift = 3;                                 % Maximum pixel drift between centroid and Gaussian fit center, default = 3
    BoxDefineLocalBkgrd = 40;                        % Box size to determine local background intensity (in pixels), default = 40
    pixel_limit_from_cy3 = 2;                        % Pixel limit between Cy5 and Cy3 molecule, default = 2 
    fold_deviation_radius = 5;                       % fold limit in radius standard deviation to exclude molecules, default = 5 
    number_Points_to_Fit = 15;                       % Number of Cy3/Cy5 pairs for LSF translation between channels, default = 15
    cy3_scale = 50000;                               % Scale for Cy3 intensity distribution to identify outliers, default = 5000
    number_of_frames_to_average = 10;                % Number of frames to average at the start of FOV to ID Cy3 molecules, default = 10
    max_distance_between_cy3_cy5 = 2;                % Maximum number of pixels between Cy3 and Cy5 molecules, default = 2

    % Create a table to store these parameters
    T = table([number_cy3_excitation_frames; number_cy5_excitation_frames; exposuret; delay;starting_frame; last_movie_frame; border; boxwidth; ...
               rectangle; HighCutoff; LowCutoff; drkwindow; area; minimum_distance_between_centers; BGcut; ...
               Gauss_Drift; BoxDefineLocalBkgrd; pixel_limit_from_cy3;fold_deviation_radius; number_Points_to_Fit; cy3_scale;number_of_frames_to_average;max_distance_between_cy3_cy5], ...
               'RowNames', {'number_cy3_excitation_frames', 'number_cy5_excitation_frames', 'exposuret', 'delay','starting_frame', 'last_movie_frame', 'border', ...
               'boxwidth', 'rectangle', 'HighCutoff', 'LowCutoff', 'drkwindow', 'area', 'minimum_distance_between_centers', ...
               'BGcut', 'Gauss_Drift', 'BoxDefineLocalBkgrd','pixel_limit_from_cy3', 'fold_deviation_radius', 'number_Points_to_Fit',...
               'cy3_scale','number_of_frames_to_average','max_distance_between_cy3_cy5'});

    % Write the table to a text file
    writetable(T, fullfile(dataout, 'parms_in.txt'), 'WriteRowNames', true, 'WriteVariableNames', false);
end
