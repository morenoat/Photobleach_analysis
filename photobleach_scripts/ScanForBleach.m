classdef ScanForBleach < handle
    properties
        UIFigure;
        Menu;
        Analysis;
        next_trajectory;
        previous_trajectory;
        remove_trajectory;
        set_y_limit;
        yLimitLabel

        UIAxesCy5;  % UIAxes for Cy5 plot
        UIAxesCy3;  % UIAxes for Cy3 plot

        UIAxesStoichiometry;  % UIAxes for bleaching distribution plot
        UIAxesInitial_state;  % UIAxes for initial state intensity distribution plot
        UIAxesFinal_state;  % UIAxes for final state intensity distribution plot
        UIAxesBleach_state;  % UIAxes for final state intensity distribution plot

        currentTrajectory = 1; % Track the current trajectory index
        parametersFilled = false; % Track if parameters are filled
        y_limit=3;%set initial limit for y axis


        % Menu items
        loadTrajectoriesMenuItem;
        SetBinValuesMenuItem;
        SaveDataMenuItem;

        % Properties to store loaded data
        molecules;
        data;
        FOV;
        dark;

        % Properties for Parameters
        parameters;
    end

    methods
        function app = ScanForBleach()
            % Define figure dimensions
            figureWidth = 1560;
            figureHeight = 790;

            % Center the figure on screen
            screenSize = get(0, 'ScreenSize');
            posX = (screenSize(3) - figureWidth) / 2;
            posY = (screenSize(4) - figureHeight) / 2;

            % Initialize the UI
            app.UIFigure = uifigure('Name', 'Scan Trajectories for Bleach Events', 'Position', [posX, posY, figureWidth, figureHeight]);
            set(app.UIFigure, 'Visible', 'on');
            app.createComponents();
        end

        function createComponents(app)
            % Add Menu and Buttons
            app.Menu = uimenu(app.UIFigure, 'Text', 'Menu');
            app.loadTrajectoriesMenuItem = uimenu(app.Menu, 'Text', 'Load Trajectories', 'Enable', 'off', 'MenuSelectedFcn', @(~,~)app.LoadTrajectoriesMenuSelected());
            uimenu(app.Menu, 'Text', 'Fill Parameters', 'MenuSelectedFcn', @(~,~)app.FillParametersMenuSelected());
            uimenu(app.Menu, 'Text', 'Load Parameters', 'MenuSelectedFcn', @(~,~)app.LoadParametersMenuSelected());
            uimenu(app.Menu, 'Text', 'Close', 'MenuSelectedFcn', @(~,~)closeApp(app));

            % Add analysis Menu
            app.Analysis = uimenu(app.UIFigure, 'Text', 'Analysis');
            uimenu(app.Analysis, 'Text', 'ID Bleaching', 'MenuSelectedFcn', @(~,~)app.BleachMenuSelected());
            app.SetBinValuesMenuItem=uimenu(app.Analysis, 'Text', 'Set Bin Parameters', 'Enable', 'off', 'MenuSelectedFcn', @(~,~)app.openBinSettingsDialog());
            % Add a new menu item for saving data in the Analysis menu
            app.SaveDataMenuItem=uimenu(app.Analysis, 'Text', 'Save Data', 'Enable', 'off','MenuSelectedFcn', @(~,~)app.SaveDataMenuSelected());


            % Add Navigation Buttons
            app.next_trajectory = uicontrol(app.UIFigure, ...
                'Style', 'pushbutton', ...
                'String', 'Next', ...
                'Position', [318, 50, 100, 30], ...
                'Callback', @(~,~)app.NextTrajectoryButtonPushed());

            app.previous_trajectory = uicontrol(app.UIFigure, ...
                'Style', 'pushbutton', ...
                'String', 'Previous', ...
                'Position', [125, 50, 100, 30], ...
                'Callback', @(~,~)app.PreviousTrajectoryButtonPushed());

            app.remove_trajectory = uicontrol(app.UIFigure, ...
                'Style', 'pushbutton', ...
                'String', 'Remove Trajectory', ...
                'Position', [536, 50, 150, 30], ...
                'Callback', @(~,~)app.RemoveTrajectoryButtonPushed());

            % Create the text label
            app.yLimitLabel = uicontrol(app.UIFigure, ...
                'Style', 'text', ...
                'String', 'Set Y limit:', ...
                'Position', [75, 100, 50, 30], ... % Adjust the position and size as needed
                'HorizontalAlignment', 'right'); % Align text to the right

            app.set_y_limit = uicontrol(app.UIFigure, ...
                'Style', 'slider', ...
                'Min', 1, ... % Set the minimum value
                'Max', 20, ... % Set the maximum value
                'Value', 3, ... % Set the initial value (optional)
                'SliderStep', [1/(20-1), 1/(20-1)], ... % Step size of 1
                'Enable', 'off', ...
                'Position', [125, 100, 150, 30], ...
                'Callback', @(src,~)app.set_y_limit_slider(src));


            % Create Axes for Plotting
            app.UIAxesCy5 = uiaxes(app.UIFigure, 'Position', [100, 150, 750, 150]);
            title(app.UIAxesCy5, 'Cy5 Intensity');
            xlabel(app.UIAxesCy5, 'Time (sec)');
            ylabel(app.UIAxesCy5, 'Intensity');

            app.UIAxesCy3 = uiaxes(app.UIFigure, 'Position', [100, 350, 750, 150]);
            title(app.UIAxesCy3, 'Cy3 Intensity');
            xlabel(app.UIAxesCy3, 'Time (sec)');
            ylabel(app.UIAxesCy3, 'Intensity');


            % Create UIAxes for initial state intensity distribution plot, initially hidden
            app.UIAxesInitial_state = uiaxes(app.UIFigure, 'Position', [100, 550, 350, 225], 'Visible', 'off');
            title(app.UIAxesInitial_state, 'Initial State Distribution');
            xlabel(app.UIAxesInitial_state, 'Intensity');
            ylabel(app.UIAxesInitial_state, 'probability density');

            % Create UIAxes for final state intensity distribution plot, initially hidden
            app.UIAxesFinal_state = uiaxes(app.UIFigure, 'Position', [500, 550, 350, 225], 'Visible', 'off');
            title(app.UIAxesFinal_state, 'Final State Distribution');
            xlabel(app.UIAxesFinal_state, 'Intensity');
            ylabel(app.UIAxesFinal_state, 'probability density');

            % Create UIAxes for bleaching state distribution plot, initially hidden
            app.UIAxesBleach_state = uiaxes(app.UIFigure, 'Position', [900, 550, 350, 225], 'Visible', 'off');
            title(app.UIAxesBleach_state, 'Bleach State Distribution');
            xlabel(app.UIAxesBleach_state, 'Intensity');
            ylabel(app.UIAxesBleach_state, 'probability density');

            % Create UIAxes for bleaching distribution plot, initially hidden
            app.UIAxesStoichiometry = uiaxes(app.UIFigure, 'Position', [1000, 150, 500, 350], 'Visible', 'off');
            title(app.UIAxesStoichiometry, 'Bleaching Distribution');
            xlabel(app.UIAxesStoichiometry, 'Stoichiometry');
            ylabel(app.UIAxesStoichiometry, 'probability density');







        end

        % Callback to load trajectories from a file
        function LoadTrajectoriesMenuSelected(app)
            % Open a file selection dialog
            [file, path] = uigetfile('*.mat', 'Select Trajectory Data File');
            if isequal(file, 0)
                disp('No file selected');
                return;
            end


            % Store data in class properties (modify as needed based on your file structure)

            app.data = load(fullfile(path, file));
            app.FOV=app.parameters.field_of_view;
            app.dark=app.parameters.dark_roi;


            if app.data.st{app.FOV}.analyzed==false
                app.data.st{app.FOV}.cy5_dark_roi_intensity_scale=app.data.st{app.FOV}.cy5_dark_roi_intensity;
                app.data.st{app.FOV}.cy3_dark_roi_intensity_scale=app.data.st{app.FOV}.cy3_dark_roi_intensity;

                app.data.st{app.FOV}.cy5_dark_roi_fit_intensity=nan(size(app.data.st{app.FOV}.cy5_dark_roi_intensity,1),size(app.data.st{app.FOV}.cy5_dark_roi_intensity,2));
                app.data.st{app.FOV}.cy3_dark_roi_intensity_smooth=nan(size(app.data.st{app.FOV}.cy5_dark_roi_intensity,1),size(app.data.st{app.FOV}.cy3_dark_roi_intensity,2));

                app.data.st{app.FOV}.cy5_intensity_scale=app.data.st{app.FOV}.cy5_intensity;
                app.data.st{app.FOV}.cy3_intensity_scale=app.data.st{app.FOV}.cy3_intensity;

                app.data.st{app.FOV}.cy5_fit_intensity=nan(size(app.data.st{app.FOV}.cy5_intensity,1),size(app.data.st{app.FOV}.cy5_intensity,2));
                app.data.st{app.FOV}.cy3_intensity_smooth=nan(size(app.data.st{app.FOV}.cy3_intensity,1),size(app.data.st{app.FOV}.cy3_intensity,2));


            end

            disp('Data loaded successfully');
            app.set_y_limit.Enable = 'on';
            app.updateTrajectoryPlot(); % Plot the first trajectory after loading
            % Bring the main GUI to the front
            app.UIFigure.Visible = 'on';  % Ensure it's visible
            uistack(app.UIFigure, 'top');  % Bring the figure to the front

        end
        % Define the callback function to clear data and close the app
        function closeApp(app)
            % Clear any data stored in the app (optional)
            app.data = [];  % Clear your data or any other properties as needed


            % Close the UIFigure
            delete(app.UIFigure);
            clc;
            % Clear the app variable from the base workspace
            evalin('base', 'clear app');
            evalin('base', 'close');
        end
        % Callback for Next Trajectory button
        function NextTrajectoryButtonPushed(app)
            %disp('Next Trajectory button pushed');
            app.currentTrajectory = app.currentTrajectory + 1;
            app.updateTrajectoryPlot();
        end

        % Callback for Previous Trajectory button
        function PreviousTrajectoryButtonPushed(app)
            %disp('Previous Trajectory button pushed');
            app.currentTrajectory = max(1, app.currentTrajectory - 1);
            app.updateTrajectoryPlot();
        end

        % Callback for Remove Trajectory button
        function RemoveTrajectoryButtonPushed(app)
            disp(['Remove Trajectory ' num2str(app.currentTrajectory)]);
            app.molecules(app.currentTrajectory,3)=0;
            if app.dark==0
                app.data.st{app.FOV}.cy5_final_spots(app.currentTrajectory,3)=0;
            else
                app.data.st{app.FOV}.cy5_dark_final_roi(app.currentTrajectory,3)=0;
            end

        end

        % Callback for set_y_limit button
        function set_y_limit_slider(app,src)
            % Update the y_limit property with the slider's current value
            app.y_limit = round(src.Value); % Round to ensure whole numbers

            % Call the plot update function with the new y-limit
            app.updateTrajectoryPlot();

        end



        % Callback to open parameter input GUI
        function FillParametersMenuSelected(app)
            parameterGui = ParameterInputGUI(app); % Pass the main GUI reference
        end

        % Callback to analyze Cy5 bleaching events
        function BleachMenuSelected(app)
            disp('Bleach Analysis selected');

            % Ensure data exists before proceeding
            if isempty(app.data.st{app.FOV})
                disp('No data loaded. Please load data first.');
                return;
            end

            [data_out]=identify_bleaching(app.data.st{app.FOV},app.parameters);
            app.data.st{app.FOV}=data_out;
            app.data.st{app.FOV}.analyzed=true;



            disp('Bleach Analysis completed');



            % Bring the main GUI to the front
            app.UIFigure.Visible = 'on';  % Ensure it's visible
            uistack(app.UIFigure, 'top');  % Bring the figure to the front
            % Alternatively, you can use:
            % app.UIFigure.WindowState = 'normal';  % Ensure it’s not minimized


        end


        % Update plot based on the current trajectory
        function updateTrajectoryPlot(app)
            n=app.FOV;
            mol_num = app.currentTrajectory;
            if isempty(app.data.st{n}.cy5_intensity) || isempty(app.data.st{n}.cy3_intensity)
                disp('No data loaded');
                return;
            end

            if app.dark==0

                app.plot_trajectory(app.UIAxesCy5, app.UIAxesCy3, ...
                    app.data.st{n}.cy5_intensity_scale, app.data.st{n}.cy5_fit_intensity, ...
                    app.data.st{n}.cy3_intensity_scale, mol_num, ...
                    app.data.st{n}.cy3_intensity_smooth,app.data.st{n}.bleach_analysis, app.data.st{n}.cy5_time, ...
                    app.data.st{n}.cy3_time, app.y_limit);
            else

                app.plot_trajectory(app.UIAxesCy5, app.UIAxesCy3, ...
                    app.data.st{n}.cy5_dark_roi_intensity_scale, app.data.st{n}.cy5_dark_roi_fit_intensity, ...
                    app.data.st{n}.cy3_dark_roi_intensity_scale, mol_num, ...
                    app.data.st{n}.cy3_dark_roi_intensity_smooth,app.data.st{n}.bleach_analysis_dark, app.data.st{n}.cy5_time, ...
                    app.data.st{n}.cy3_time, app.y_limit);

            end
        end

        function plot_trajectory(app, ax1, ax2, cy5_intensity, cy5_fit_intensity, cy3_intensity, mol_num, cy3_smoothed_intensity,bleach_analysis, cy5_time, cy3_time, y_limit)
            % Clear existing plots
            cla(ax1);
            cla(ax2);

            if app.data.st{app.FOV}.analyzed==true
                % Plot Cy5 Data
                plot(ax1, cy5_time, cy5_intensity(:, mol_num), 'r', 'LineWidth', 3); hold(ax1, 'on');
                plot(ax1, cy5_time, cy5_fit_intensity(:, mol_num), 'k', 'LineWidth', 2);
                ax1.YLim = [-0.2 y_limit];
                xRound = ceil(cy5_time(end,1)/10)*10;
                ax1.XLim = [0 xRound];
                ax1.XTick = 0:(ceil((xRound/5)/10)*10):xRound;
                if ~isempty(bleach_analysis{mol_num}.changepoint_list)

                    bleaching_points(:,1)=cy5_time(bleach_analysis{mol_num}.changepoint_list(:,2),1);
                    bleaching_points(:,2)=cy5_fit_intensity(bleach_analysis{mol_num}.changepoint_list(:,2),mol_num);
                else
                    bleaching_points(:,1)=0;
                    bleaching_points(:,2)=0;
                end
                scatter(ax1,bleaching_points(:,1),bleaching_points(:,2),100,'b','filled');

                title(ax1, sprintf('Cy5 - Molecule %d', mol_num));

                % Plot Cy3 Data
                plot(ax2, cy3_time, cy3_intensity(:, mol_num), 'g', 'LineWidth', 3); hold(ax2, 'on');
                plot(ax2, cy3_time, cy3_smoothed_intensity(:, mol_num), 'k', 'LineWidth', 2);
                ax2.YLim = [-0.2 y_limit];
                xRound = ceil(cy3_time(end,1)/10)*10;
                ax2.XLim = [0 xRound];
                ax2.XTick = 0:(ceil((xRound/5)/10)*10):xRound;


                title(ax2, sprintf('Cy3 - Molecule %d', mol_num));



                % Make the UIAxes for histograms and bleaching distribution visible
                app.UIAxesInitial_state.Visible = 'on';
                app.UIAxesFinal_state.Visible = 'on';
                app.UIAxesBleach_state.Visible = 'on';
                app.UIAxesStoichiometry.Visible = 'on';


                cla(app.UIAxesInitial_state);
                cla(app.UIAxesFinal_state);
                cla(app.UIAxesBleach_state);
                cla(app.UIAxesStoichiometry);

                if app.dark==0

                    bleach_data=app.data.st{app.FOV}.cy5_tabulate_bleaching.list_bleaching_events;
                else
                    bleach_data=app.data.st{app.FOV}.cy5_dark_roi_tabulate_bleaching.list_bleaching_events;

                end


                % Plot histograms
                % Define the bin edges
                binEdgesInitial = 0:0.5:5;

                % Plot the histogram for initial state
                histogram(app.UIAxesInitial_state, bleach_data(:,2), binEdgesInitial, 'Normalization', 'pdf');

                % Plot the histogram for final state
                histogram(app.UIAxesFinal_state, bleach_data(:,3), binEdgesInitial, 'Normalization', 'pdf');

                % Plot the histogram for bleach state
                histogram(app.UIAxesBleach_state, bleach_data(:,4), binEdgesInitial, 'Normalization', 'pdf');

                % Define the bin edges for stoichiometry
                binEdgesStoichiometry = 0:1:10;

                % Plot the histogram for stoichiometry
                histogram(app.UIAxesStoichiometry, bleach_data(:,1), binEdgesStoichiometry, 'Normalization', 'pdf');

                app.SaveDataMenuItem.Enable = 'on';
                app.SetBinValuesMenuItem.Enable = 'on';







            else
                cy5_norm=median(max(cy5_intensity,[],1));
                cy3_norm=median(max(cy3_intensity,[],1));

                tmp_cy5=(cy5_intensity(:, mol_num)./cy5_norm);

                if prctile(tmp_cy5,95)>app.y_limit
                    tmp_cy5=tmp_cy5./prctile(tmp_cy5,95);
                end

                % Plot Cy5 Data
                plot(ax1, cy5_time,tmp_cy5, 'r', 'LineWidth', 3); hold(ax1, 'on');
                ax1.YLim = [-0.2 y_limit];
                xRound = ceil(cy5_time(end,1)/10)*10;
                ax1.XLim = [0 xRound];
                ax1.XTick = 0:(ceil((xRound/5)/10)*10):xRound;


                title(ax1, sprintf('Cy5 - Molecule %d', mol_num));

                tmp_cy3=(cy3_intensity(:, mol_num)./cy3_norm);

                if prctile(tmp_cy3,95)>app.y_limit
                    tmp_cy3=tmp_cy3./prctile(tmp_cy3,95);
                end

                % Plot Cy3 Data
                plot(ax2, cy3_time, tmp_cy3, 'g', 'LineWidth', 2); hold(ax2, 'on');
                ax2.YLim = [-0.2 y_limit];
                xRound = ceil(cy3_time(end,1)/10)*10;
                ax2.XLim = [0 xRound];
                ax2.XTick = 0:(ceil((xRound/5)/10)*10):xRound;


                title(ax2, sprintf('Cy3 - Molecule %d', mol_num));

            end
        end
        % Callback function for saving data
        function SaveDataMenuSelected(app)
            % Check if data exists to be saved
            if isempty(app.data)
                disp('No data to save. Please complete the analysis first.');
                return;
            end

            % Open a save dialog and save the data
            st=app.data.st;
            uisave('st', 'BleachAnalysisData');
            disp('Data saved successfully.');

            % Extract the matrix from your data structure
            if app.dark==0
                matrixData = app.data.st{app.FOV}.cy5_tabulate_bleaching.list_bleaching_events;

                % Convert the matrix to a table and assign column labels
                dataTable = array2table(matrixData, 'VariableNames', {'Stoichiometry', 'InitialIntensity', 'FinalIntensity', 'BleachIntensity'});

                % Specify the file name for the CSV file
                fileName = 'bleaching_cy5_data.csv';
            else

                matrixData = app.data.st{app.FOV}.cy5_dark_roi_tabulate_bleaching.list_bleaching_events;

                % Convert the matrix to a table and assign column labels
                dataTable = array2table(matrixData, 'VariableNames', {'Stoichiometry', 'InitialIntensity', 'FinalIntensity', 'BleachIntensity'});

                % Specify the file name for the CSV file
                fileName = 'bleaching_cy5_dark_roi_data.csv';



            end

            % Save the table as a CSV file
            writetable(dataTable, fileName);

            % Notify user (optional)
            disp(['Data saved as ', fileName]);
        end






    end

    % Method for enabling load trajectories menu after parameters are filled
    methods
        function enableLoadTrajectories(app)
            app.parametersFilled = true;
            app.loadTrajectoriesMenuItem.Enable = 'on';
        end
    end

    methods
        % Method to load parameters from a .txt file
        function LoadParametersMenuSelected(app)
            % Open a file selection dialog
            % [file, path] = uigetfile('*.txt', 'Select Parameters File');
            % Bring the main application window to the front
            drawnow;  % Ensure all pending GUI operations are completed
            figure(app.UIFigure);  % Bring the main figure to the front

            % Call uigetfile
            [file, path] = uigetfile('*.txt', 'Select Parameters File');
            if isequal(file, 0)
                disp('No file selected');
                return;
            end

            % Read the table from the file
            paramTable = readtable(fullfile(path, file), 'Delimiter', '\t');


            % Iterate through the table and assign the values to properties
            for i = 1:height(paramTable)
                paramName = paramTable.ParameterName{i};
                paramValue = paramTable.Value(i);
                app.parameters=setfield(app.parameters,paramName,paramValue);

                %                 % Check if the property exists in the app and update its value
                %                 if isprop(app, paramName)
                %                     app.(paramName) = paramValue;
                %                 end
            end

            disp('Parameters loaded successfully');

            app.parametersFilled = true;
            app.loadTrajectoriesMenuItem.Enable = 'on';
            % Bring the main GUI to the front
            app.UIFigure.Visible = 'on';  % Ensure it's visible
            uistack(app.UIFigure, 'top');  % Bring the figure to the front
            % Alternatively, you can use:
            % app.UIFigure.WindowState = 'normal';  % Ensure it’s not minimized

        end



        % Method to create bin parameter dialog
        function openBinSettingsDialog(app)


            figWidth = 600;
            figHeight = 300;
            fieldHeight = 22;
            marginTop = 50;
            marginLeft = 20;
            columnSpacing = (figWidth - 2 * marginLeft) / 2;
            % Center the figure on screen
            screenSize = get(0, 'ScreenSize');
            posX = (screenSize(3) - figWidth) / 2;
            posY = (screenSize(4) - figHeight) / 2;

            % Initialize the UI figure
            d = uifigure('Name', 'Histogram Bin Parameters', 'Position', [posX, posY, figWidth, figHeight]);

            % Bin labels and default values
            BinLabels = {'Stoichiometry Bin Width:', 'Intensity Bin Width:','Stoichiometry Min Value:',...
                'Intensity Min Value:','Stoichiometry Max Value:','Intensity Max Value:'};

            defaultValues = [1, 0.5, 0, 0, 5, 5];

            % Create histogram bin fields in two columns


            for i = 1:length(BinLabels)
                row = ceil(i / 2); % Calculate row for each parameter
                col = mod(i - 1, 2); % 0 for left column, 1 for right column
                xPos = marginLeft + col * columnSpacing;
                yPos = figHeight - marginTop - row * (fieldHeight + 10); % Dynamic y position

                % Create label
                labelWidth = 200; % Width of the label
                fieldOffset = 150; % Increase this value for more space between label and field
                uilabel(d, 'Text', BinLabels{i}, ...
                    'Position', [xPos, yPos, labelWidth, fieldHeight]);

                % Create corresponding input field
                if islogical(defaultValues(i))  % Checkbox for logical (boolean) values
                    field = uicheckbox(d, 'Position', [xPos + fieldOffset, yPos, 50, fieldHeight], ...
                        'Value', defaultValues(i));
                else  % Numeric field for other values
                    field = uieditfield(d, 'numeric', ...
                        'Position', [xPos + fieldOffset, yPos, 100, fieldHeight], ...
                        'Value', defaultValues(i));
                end

                BinLabels{i}=field;


            end


            uibutton(d, 'Text', 'Apply', 'Position', [250, 25, 100, 30], 'ButtonPushedFcn', @(~, ~) applyBinSettings(BinLabels));



            function applyBinSettings(Bin_setting)

                % Define the new bin edges based on the user input
                intensity_binEdges = Bin_setting{4}.Value:Bin_setting{2}.Value: Bin_setting{6}.Value;


                % Define the new bin edges based on the user input
                stochiometry_binEdges = Bin_setting{3}.Value:Bin_setting{1}.Value: Bin_setting{5}.Value;


                % Display the settings (optional)
                disp([' Intensity Bin width set to: ', num2str(Bin_setting{2}.Value)]);
                disp(['Range: [', num2str(Bin_setting{4}.Value), ', ', num2str(Bin_setting{6}.Value), ']']);

                disp([' Stoichiometry Bin width set to: ', num2str(Bin_setting{1}.Value)]);
                disp(['Range: [', num2str(Bin_setting{3}.Value), ', ', num2str(Bin_setting{5}.Value), ']']);

                % Update the histograms with the new bin edges
                cla(app.UIAxesInitial_state);
                cla(app.UIAxesFinal_state);
                cla(app.UIAxesBleach_state);
                cla(app.UIAxesStoichiometry);

                if app.dark==0

                    bleach_data = app.data.st{app.FOV}.cy5_tabulate_bleaching.list_bleaching_events;

                else
                    bleach_data = app.data.st{app.FOV}.cy5_dark_roi_tabulate_bleaching.list_bleaching_events;
                end

                % Plot updated histograms
                histogram(app.UIAxesInitial_state, bleach_data(:,2), intensity_binEdges, 'Normalization', 'pdf');
                histogram(app.UIAxesFinal_state, bleach_data(:,3), intensity_binEdges, 'Normalization', 'pdf');
                histogram(app.UIAxesBleach_state, bleach_data(:,4), intensity_binEdges, 'Normalization', 'pdf');
                histogram(app.UIAxesStoichiometry, bleach_data(:,1), stochiometry_binEdges, 'Normalization', 'pdf');


                % Close the dialog box
                close(d);
            end

        end
    end



end







