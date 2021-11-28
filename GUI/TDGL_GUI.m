classdef TDGL_GUI < GUI_B_boundary
    properties 
        
        
    end
    
    methods
        function app = TDGL_GUI()
            app = app@GUI_B_boundary();
            
        end
    end
    
    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            trapezoidal_test_GUI(app);
        end

        % Window button down function: UIFigure
        function UIFigureWindowButtonDown(app, event)
            P = get(app.UIAxes,'CurrentPoint');
            disp(P)
        end
    end

    
end
            