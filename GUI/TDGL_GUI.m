classdef TDGL_GUI < GUI_B_boundary
    properties 
        % Startup values
        Nx
        Ny
        Nz
        Bz
        BzButtonPushed
        kappaIO
        stop
    end
    
    methods
        function app = TDGL_GUI()
            app = app@GUI_B_boundary();
            
        end
        
            
    end
    
    
    
    
%% COPIED FROM APP DESIGNER
  % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            app.StartButton.Enable = 'off';
            app.stop = 0;
            app.Nx = app.NxSpinner.Value;
            app.Ny = app.NySpinner.Value;
            app.Nz = app.NzSpinner.Value;
            app.Bz = 0;
            if app.SCTypeButtonGroup.Buttons(1).Value == 1
                app.kappaIO = 5;
            else
                app.kappaIO = 0.5;
            end
            
            app.StatusEditField.Value = "Simulation Started";
            trapezoidal_test_GUI(app);

            
            
        end

        % Window button down function: UIFigure
        function UIFigureWindowButtonDown(app, event)
            P = get(app.UIAxes,'CurrentPoint');

        end

        % Button pushed function: ApplyFIeldButton
        function ApplyFIeldButtonPushed(app, event)
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

        end

        % Value changed function: NxSpinner
        function NxSpinnerValueChanged(app, event)
            value = app.NxSpinner.Value;
            app.NySpinner.Value = value;
        end
    end



    
%%    
end
            