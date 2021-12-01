classdef TDGL_GUI < GUI_B_boundary
    properties 
        % Startup values
        Bz
        BzClick
        BzButtonPushed
        cord
        click
        kappaIO
        Nx
        Ny
        Nz
        stop
    end
    

    methods
        function app = TDGL_GUI()
            app = app@GUI_B_boundary();
        end
    end
    
    
    
    
%% COPIED FROM APP DESIGNER
    % Callbacks that handle component events
    methods 

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
 
            app.StartButton.Enable = 'off';
            cla(app.UIAxes)
            cla(app.UIAxes2)
            app.stop = 0;
            app.Bz = 0;

            app.Nx = app.NxSpinner.Value;
            app.Ny = app.NySpinner.Value;
            app.Nz = app.NzSpinner.Value;
            
            if app.SCTypeButtonGroup.Buttons(1).Value == 1
                app.kappaIO = 5;
            else
                app.kappaIO = 0.5;
            end
            
            app.StatusEditField.Value = "Simulation Started";

            trapezoidal_test_GUI(app, event);

        end

        % Button pushed function: ApplyFIeldEverywhereButton
        function ApplyFIeldEverywhereButtonPushed(app, event)
            app.Bz = app.BzSlider.Value;
            app.BzButtonPushed = 1;
        end

        % Button pushed function: ResetFieldButton
        function ResetFieldButtonPushed(app, event)
            app.BzSlider.Value = 0;
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
            app.stop = 1;
            app.StatusEditField.Value = "Simulation Stopped";
            app.StartButton.Enable = 'on';
        end

        % Value changed function: NxSpinner
        function NxSpinnerValueChanged(app, event)
            value = app.NxSpinner.Value;
            app.NySpinner.Value = value;
        end

        % Button down function: UIFigure
        function UIFigureButtonDown(app, event)
            if ~isempty(app.UIAxes.CurrentPoint)
                app.cord = get(app.UIAxes,'CurrentPoint');
            end
            app.click = 1;  
        end
    end


%%    
end
            