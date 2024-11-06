classdef ParameterInputGUI < handle
    properties
        UIFigure;
        mainApp; % Reference to main GUI



        % Parameters to be set in GUI
        scale_valueField;
        smoothingField;
        minimum_intensityField;
        fractional_changeField;
        relative_to_max_intensityField;
        move_median_sampleField;
        minimum_pointsField;
        replace_negativeField;
        dark_roiField;
        field_of_viewField;



    end

    methods
        function app = ParameterInputGUI(mainApp)
            app.mainApp = mainApp;
            app.createComponents();
        end

        function createComponents(app)
            % Set figure dimensions and spacing settings
            figWidth = 1020;
            figHeight = 400;
            fieldHeight = 22;
            marginTop = 20;
            marginLeft = 20;
            columnSpacing = (figWidth - 2 * marginLeft) / 2;
            % Center the figure on screen
            screenSize = get(0, 'ScreenSize');
            posX = (screenSize(3) - figWidth) / 2;
            posY = (screenSize(4) - figHeight) / 2;

            % Initialize the UI figure
            app.UIFigure = uifigure('Name', 'Fill Parameters', 'Position', [posX, posY, figWidth, figHeight]);

            % Parameter labels and default values
            paramNames = {'scale_value', 'smoothing', 'minimum_intensity', ...
                'fractional_change', 'relative_to_max_intensity', ...
                'move_median_sample', 'minimum_points', 'replace_negative','dark_roi', 'field_of_view'};

            defaultValues = [3, 5, 0.2, 0.35, 2, 2, 2, 0.03, false, 1];

            % Create parameter fields in two columns



            for i = 1:length(paramNames)
                row = ceil(i / 2); % Calculate row for each parameter
                col = mod(i - 1, 2); % 0 for left column, 1 for right column
                xPos = marginLeft + col * columnSpacing;
                yPos = figHeight - marginTop - row * (fieldHeight + 10); % Dynamic y position

                % Create label
                labelWidth = 200; % Width of the label
                fieldOffset = 200; % Increase this value for more space between label and field
                uilabel(app.UIFigure, 'Text', paramNames{i}, ...
                    'Position', [xPos, yPos, labelWidth, fieldHeight]);

                % Create corresponding input field
                if islogical(defaultValues(i))  % Checkbox for logical (boolean) values
                    field = uicheckbox(app.UIFigure, 'Position', [xPos + fieldOffset, yPos, 50, fieldHeight], ...
                        'Value', defaultValues(i));
                else  % Numeric field for other values
                    field = uieditfield(app.UIFigure, 'numeric', ...
                        'Position', [xPos + fieldOffset, yPos, 100, fieldHeight], ...
                        'Value', defaultValues(i));
                end

                % Assign the field to the corresponding property
                app.(paramNames{i} + "Field") = field;
            end



            % Add Submit Button
            uibutton(app.UIFigure, 'Text', 'Submit', ...
                'Position', [figWidth / 2 - 50, marginTop / 2, 100, fieldHeight], ...
                'ButtonPushedFcn', @(~, ~) app.submitParameters());
        end

        function submitParameters(app)
            % Assign parameters to main app properties
       
            parameters.scale_value = app.scale_valueField.Value;
            parameters.smoothing = app.smoothingField.Value;
            parameters.minimum_intensity = app.minimum_intensityField.Value;
            parameters.fractional_change = app.fractional_changeField.Value;
            parameters.relative_to_max_intensity = app.relative_to_max_intensityField.Value;
            parameters.move_median_sample = app.move_median_sampleField.Value;
            parameters.minimum_points = app.minimum_pointsField.Value;
            parameters.replace_negative = app.replace_negativeField.Value;
            parameters.dark_roi = app.dark_roiField.Value;
            parameters.field_of_view = app.field_of_viewField.Value;
            app.mainApp.parameters=parameters;


            %save paramters for loading at a later date

            app.saveParametersToFile()

            % Enable load trajectories menu item
            app.mainApp.enableLoadTrajectories();

            % Close the parameter input GUI
            delete(app.UIFigure);
        end
    end
    methods
        % Method to save parameters to a .txt file
        function saveParametersToFile(app)
            selpath = uigetdir ;
            filename='parameters.txt';
            % Get parameter names and values
            paramNames = fieldnames(app);
            paramValues = [];

            % Loop through each field to collect the value if it's an input field
            for i = 1:length(paramNames)
                if endsWith(paramNames{i}, 'Field')  % Check if it's an input field property
                    value = app.(paramNames{i}).Value;
                    newStr = erase( paramNames{i} , 'Field' );
                    %paramValues = [paramValues; {paramNames{i}, value}];  % Store name and value
                    paramValues = [paramValues; {newStr, value}];  % Store name and value
                end
            end

            % Convert to table
            paramTable = cell2table(paramValues, 'VariableNames', {'ParameterName', 'Value'});

            % Write to .txt file
            writetable(paramTable, fullfile(selpath,filename), 'Delimiter', '\t');

            disp(['Parameters saved to ', filename]);
        end
    end


end
