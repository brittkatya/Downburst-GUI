classdef OBV_AbdElaal_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        Panel_3                         matlab.ui.container.Panel
        FixedParametersCheckBox         matlab.ui.control.CheckBox
        PlottingOptionsPanel            matlab.ui.container.Panel
        PlotverticalvelocityCheckBox    matlab.ui.control.CheckBox
        PlothorizontalvelocityCheckBox  matlab.ui.control.CheckBox
        PlotmagnitudeCheckBox           matlab.ui.control.CheckBox
        Panel                           matlab.ui.container.Panel
        zField                          matlab.ui.control.NumericEditField
        zSlider                         matlab.ui.control.Slider
        zmLabel                         matlab.ui.control.Label
        Panel_2                         matlab.ui.container.Panel
        rSlider                         matlab.ui.control.Slider
        rField                          matlab.ui.control.NumericEditField
        zmLabel_2                       matlab.ui.control.Label
        ExportDataTableButton_2         matlab.ui.control.Button
        ExportPlotButton_2              matlab.ui.control.Button
        VerticalProfileModelParametersPanel  matlab.ui.container.Panel
        MaxhorizontalvelField_2         matlab.ui.control.NumericEditField
        MaxhorizontalvelSlider_2        matlab.ui.control.Slider
        MaximumhorizontalvelocitymsLabel_2  matlab.ui.control.Label
        JetDiameterField_2              matlab.ui.control.NumericEditField
        JetDiameterSlider_2             matlab.ui.control.Slider
        JetDiameterkmLabel_2            matlab.ui.control.Label
        MaxverticalvelField             matlab.ui.control.NumericEditField
        MaxverticalvelSlider            matlab.ui.control.Slider
        MaximumverticalvelocitymsLabel  matlab.ui.control.Label
        ExportDataTableButton           matlab.ui.control.Button
        ExportPlotButton                matlab.ui.control.Button
        RadialProfileModelParametersPanel  matlab.ui.container.Panel
        MaxverticalvelField_2           matlab.ui.control.NumericEditField
        MaxverticalvelSlider_2          matlab.ui.control.Slider
        MaximumverticalvelocitymsLabel_2  matlab.ui.control.Label
        JetDiameterField                matlab.ui.control.NumericEditField
        JetDiameterSlider               matlab.ui.control.Slider
        JetDiameterkmLabel              matlab.ui.control.Label
        MaxhorizontalvelField           matlab.ui.control.NumericEditField
        MaxhorizontalvelSlider          matlab.ui.control.Slider
        MaximumhorizontalvelocitymsLabel  matlab.ui.control.Label
        AbdElaalMillsMa2013Label        matlab.ui.control.Label
        VerticalAxes                    matlab.ui.control.UIAxes
        HorizontalAxes                  matlab.ui.control.UIAxes
    end

    
    methods (Access = private)
        
        function u_r = compute_ur(app, r, sf_u, ur_1, ur_2)
            % Horizontal downburst velocity for radial profile
            u_r = (sf_u .* r)./2 .* ur_1 .* ur_2;
        end
        
        function w_r = compute_wr(app, sf_w, wr_1, wr_2)
            % Vertical downburst velocity for radial profile
            w_r = -sf_w .* wr_1 .* wr_2;
        end

        function u = compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2)
            u_r = (sf_u .* r)./2 .* ur_1 .* ur_2;
            w_r = -sf_w .* wr_1 .* wr_2;

            % Downburst velocity magnitude for radial profile
            u = sqrt(u_r.^2 + w_r.^2);
        end

        function u_z = compute_uz(app, r, sf_u, uz_1, uz_2)
            % Horizontal downburst velocity for radial profile
            u_z = (sf_u .* r)./2 .* uz_1 .* uz_2;
        end

        function w_z = compute_wz(app, sf_w, wz_1, wz_2)
            % Vertical downburst velocity for vertical profile
            w_z = -sf_w .* wz_1 .* wz_2;
        end

        function w = compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2)
            u_z = (sf_u .* r)./2 .* uz_1 .* uz_2;
            w_z = -sf_w .* wz_1 .* wz_2;

            % Downburst velocity magnitude for radial profile
            w = sqrt(u_z.^2 + w_z.^2);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Set up plotting options menu
            app.PlotmagnitudeCheckBox.Value = 1;
            app.PlothorizontalvelocityCheckBox.Value = 1;
            app.PlotverticalvelocityCheckBox.Value = 1;
            app.FixedParametersCheckBox.Value = 1;            

            % Define parameters for radial profile
            r = 0:0.1:10000;
            u_max = app.MaxhorizontalvelField.Value;
            w_max = app.MaxverticalvelField_2.Value;
            D = app.JetDiameterField.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            z = app.zField.Value;
%             gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
%             de = 2.1;
%             eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
%             kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
%             ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%             c1 = -0.16;
%             c2 = 1.19;
            de = 2.1;
            gam = 0.35;
            eps = 0.18;
            kap = 0.06;
            ch = 1;
            c1 = -4.1;
            c2 = 1.1;
            ps_r = (de*r.^2./r_max^2).^gam;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            ur_1 = exp(-(2*gam - (de.*r.^2./r_max^2).^gam).^2) + eps .* exp(-kap.*(r.^2./r_max^2));
            ur_2 = (z/D)^0.1 * exp(c1*(z/D)^1.1);
            wr_1 = (1 + 2.*gam.*ps_r.*(2*gam-ps_r)) .* exp(-(2*gam - ps_r).^2) + eps .* exp(-kap.*(r.^2./r_max^2).^ch) .* (1 - kap.*ch.*(r./r_max).^(2*ch));
            wr_2 = (D/-4.51) * (exp(-4.51*(z/D).^1.1)-1);

            % Plot radial profile horizontal and vertical velocities as well as magnitude
            plot(app.HorizontalAxes, r, compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2), 'k', r, compute_ur(app, r, sf_u, ur_1, ur_2), 'r', r, compute_wr(app, sf_w, wr_1, wr_2), 'b');
            p = app.HorizontalAxes.Children;
            p(1).LineWidth = 1.5;
            p(2).LineWidth = 1.5;
            p(3).LineWidth = 1.5;
            legend(app.HorizontalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
            legend(app.HorizontalAxes, 'Location', 'southoutside');            

%             % Adjust axis for radial profile
%             r_lim = 2 * r_max;
%             xlim(app.HorizontalAxes, [0, (r_lim+2000)]);

            % Define parameters for vertical profile
            z = 0:0.1:1000;
            u_max = app.MaxhorizontalvelField_2.Value;
            w_max = app.MaxverticalvelField.Value;
            D = app.JetDiameterField_2.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            r = app.rField.Value * 1000;
            ps_z = (2*r.^2./D^2).^0.35;
            % These parameters are from equations 29-32, but they make the
            % vertical profile look odd
            gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
            de = 2.1;
            eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
            kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
            ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%             % These parameters are from paragraph 2 of section 4.4 and makes the vertical plot look better but
%             % they're for horizontal profile
%             de = 2.1;
%             gam = 0.35;
%             eps = 0.18;
%             kap = 0.06;
%             ch = 1;
            c1 = -0.16;
            c2 = 1.19;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            uz_1 = exp(-(0.7 - (2.*r.^2./D^2).^0.35).^2) + 0.18*exp(-0.06.*(r./D).^2);
            uz_2 = (z./z_max).^(c2-1) .* exp(c1.*(z./z_max).^c2);
            wz_1 = (1 + 0.7*ps_z*(0.7-ps_z)) * exp(-(0.7 - ps_z)^2) + 0.18*exp(-0.06*(r/D)^2) .* (1 - kap.*ch.*((r/D)^2));
            wz_2 = (z_max/(c1*c2)) * (exp(c1.*(z./z_max).^c2)-1);

            % Plot vertical profile horizontal and vertical velocities as well as magnitude
            plot(app.VerticalAxes, compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2), z, 'k', compute_uz(app, r, sf_u, uz_1, uz_2), z, 'r', compute_wz(app, sf_w, wz_1, wz_2), z, 'b');
            p = app.VerticalAxes.Children;
            p(1).LineWidth = 1.5;
            p(2).LineWidth = 1.5;
            p(3).LineWidth = 1.5;
            legend(app.VerticalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
            legend(app.VerticalAxes, 'Location', 'southoutside');
        end

        % Value changed function: JetDiameterField, MaxhorizontalvelField, 
        % MaxverticalvelField_2, zField
        function FieldValueChanged(app, event)
            % Define parameters for radial profile
            r = 0:0.1:10000;
            u_max = app.MaxhorizontalvelField.Value;
            w_max = app.MaxverticalvelField_2.Value;
            D = app.JetDiameterField.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            z = app.zField.Value;
            de = 2.1;
            gam = 0.35;
            eps = 0.18;
            kap = 0.06;
            ch = 1;
            c1 = -4.1;
            c2 = 1.1;
            ps_r = (de*r.^2./r_max^2).^gam;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            ur_1 = exp(-(2*gam - (de.*r.^2./r_max^2).^gam).^2) + eps .* exp(-kap.*(r.^2./r_max^2));
            ur_2 = (z/D)^0.1 * exp(c1*(z/D)^1.1);
            wr_1 = (1 + 2.*gam.*ps_r.*(2*gam-ps_r)) .* exp(-(2*gam - ps_r).^2) + eps .* exp(-kap.*(r.^2./r_max^2).^ch) .* (1 - kap.*ch.*(r./r_max).^(2*ch));
            wr_2 = (D/-4.51) * (exp(-4.51*(z/D).^1.1)-1);

            % Plot radial profile horizontal and vertical velocities as well as magnitude
            plot(app.HorizontalAxes, r, compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2), 'k', r, compute_ur(app, r, sf_u, ur_1, ur_2), 'r', r, compute_wr(app, sf_w, wr_1, wr_2), 'b');
            p = app.HorizontalAxes.Children;
            p(1).LineWidth = 1.5;
            p(2).LineWidth = 1.5;
            p(3).LineWidth = 1.5;
            legend(app.HorizontalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
            legend(app.HorizontalAxes, 'Location', 'southoutside');            

%             % Adjust axis for radial profile
%             r_lim = 2 * r_max;
%             xlim(app.HorizontalAxes, [0, (r_lim+2000)]);

            % "Plotting options" menu check box adjustments
            if app.PlotmagnitudeCheckBox.Value == 1
                p(3).Visible = 'on';
            else p(3).Visible = 'off';
            end
            if app.PlothorizontalvelocityCheckBox.Value == 1
                p(2).Visible = 'on';
            else p(2).Visible = 'off';
            end
            if app.PlotverticalvelocityCheckBox.Value == 1
                p(1).Visible = 'on';
            else p(1).Visible = 'off';
            end

            % Update sliders to match field values
            if app.MaxhorizontalvelField.Value ~= app.MaxhorizontalvelSlider.Value
                app.MaxhorizontalvelSlider.Value = app.MaxhorizontalvelField.Value;
            else 
            end
            if app.MaxverticalvelField_2.Value ~= app.MaxverticalvelSlider_2.Value
                app.MaxverticalvelSlider_2.Value = app.MaxverticalvelField_2.Value;
            end
            if app.JetDiameterField.Value ~= app.JetDiameterSlider.Value
                app.JetDiameterSlider.Value = app.JetDiameterField.Value;
            else
            end
            if app.zField.Value ~= app.zSlider.Value
                app.zSlider.Value = app.zField.Value;
            else 
            end

            % Update field value if check box is checked
            if app.FixedParametersCheckBox.Value == 1
                app.MaxhorizontalvelField_2.Value = app.MaxhorizontalvelField.Value;
                app.MaxhorizontalvelSlider_2.Value = app.MaxhorizontalvelSlider.Value;
                app.MaxverticalvelField.Value = app.MaxverticalvelField_2.Value;
                app.MaxverticalvelSlider.Value = app.MaxverticalvelSlider_2.Value;
                app.JetDiameterField_2.Value = app.JetDiameterField.Value;
                app.JetDiameterSlider_2.Value = app.JetDiameterSlider.Value;

                % Define parameters for vertical profile
                z = 0:0.1:1000;
                u_max = app.MaxhorizontalvelField_2.Value;
                w_max = app.MaxverticalvelField.Value;
                D = app.JetDiameterField_2.Value * 1000;
                r_max = 1.125*D;
                z_max = 0.025*D;
                r = app.rField.Value * 1000;
                ps_z = (2*r.^2./D^2).^0.35;
                gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
                de = 2.1;
                eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
                kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
                ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%                 de = 2.1;
%                 gam = 0.35;
%                 eps = 0.18;
%                 kap = 0.06;
%                 ch = 1;
                c1 = -0.16;
                c2 = 1.19;
                sf_u = u_max / (0.338 * r_max);
                sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));
    
                % Individual parts of the functions
                uz_1 = exp(-(0.7 - (2.*r.^2./D^2).^0.35).^2) + 0.18*exp(-0.06.*(r./D).^2);
                uz_2 = (z./z_max).^(c2-1) .* exp(c1.*(z./z_max).^c2);
                wz_1 = (1 + 0.7*ps_z*(0.7-ps_z)) * exp(-(0.7 - ps_z)^2) + 0.18*exp(-0.06*(r/D)^2) .* (1 - kap.*ch.*((r/D)^2));
                wz_2 = (z_max/(c1*c2)) * (exp(c1.*(z./z_max).^c2)-1);
    
                % Plot vertical profile horizontal and vertical velocities as well as magnitude
                plot(app.VerticalAxes, compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2), z, 'k', compute_uz(app, r, sf_u, uz_1, uz_2), z, 'r', compute_wz(app, sf_w, wz_1, wz_2), z, 'b');
                p = app.VerticalAxes.Children;
                p(1).LineWidth = 1.5;
                p(2).LineWidth = 1.5;
                p(3).LineWidth = 1.5;
                legend(app.VerticalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
                legend(app.VerticalAxes, 'Location', 'southoutside');

                % "Plotting options" menu check box adjustments
                if app.PlotmagnitudeCheckBox.Value == 1
                    p(3).Visible = 'on';
                else p(3).Visible = 'off';
                end
                if app.PlothorizontalvelocityCheckBox.Value == 1
                    p(2).Visible = 'on';
                else p(2).Visible = 'off';
                end
                if app.PlotverticalvelocityCheckBox.Value == 1
                    p(1).Visible = 'on';
                else p(1).Visible = 'off';
                end
            else
            end
        end

        % Value changed function: JetDiameterSlider, 
        % MaxhorizontalvelSlider, MaxverticalvelSlider_2, zSlider
        function SliderValueChanged(app, event)
            % Define parameters for radial profile
            r = 0:0.1:10000;
            u_max = app.MaxhorizontalvelSlider.Value;
            w_max = app.MaxverticalvelSlider_2.Value;
            D = app.JetDiameterSlider.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            z = app.zSlider.Value;
            de = 2.1;
            gam = 0.35;
            eps = 0.18;
            kap = 0.06;
            ch = 1;
            c1 = -4.1;
            c2 = 1.1;
            ps_r = (de*r.^2./r_max^2).^gam;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            ur_1 = exp(-(2*gam - (de.*r.^2./r_max^2).^gam).^2) + eps .* exp(-kap.*(r.^2./r_max^2));
            ur_2 = (z/D)^0.1 * exp(c1*(z/D)^1.1);
            wr_1 = (1 + 2.*gam.*ps_r.*(2*gam-ps_r)) .* exp(-(2*gam - ps_r).^2) + eps .* exp(-kap.*(r.^2./r_max^2).^ch) .* (1 - kap.*ch.*(r./r_max).^(2*ch));
            wr_2 = (D/-4.51) * (exp(-4.51*(z/D).^1.1)-1);

            % Plot radial profile horizontal and vertical velocities as well as magnitude
            plot(app.HorizontalAxes, r, compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2), 'k', r, compute_ur(app, r, sf_u, ur_1, ur_2), 'r', r, compute_wr(app, sf_w, wr_1, wr_2), 'b');
            p = app.HorizontalAxes.Children;
            p(1).LineWidth = 1.5;
            p(2).LineWidth = 1.5;
            p(3).LineWidth = 1.5;
            legend(app.HorizontalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
            legend(app.HorizontalAxes, 'Location', 'southoutside');            

%             % Adjust axis for radial profile
%             r_lim = 2 * r_max;
%             xlim(app.HorizontalAxes, [0, (r_lim+2000)]);

            % "Plotting options" menu check box adjustments
            if app.PlotmagnitudeCheckBox.Value == 1
                p(3).Visible = 'on';
            else p(3).Visible = 'off';
            end
            if app.PlothorizontalvelocityCheckBox.Value == 1
                p(2).Visible = 'on';
            else p(2).Visible = 'off';
            end
            if app.PlotverticalvelocityCheckBox.Value == 1
                p(1).Visible = 'on';
            else p(1).Visible = 'off';
            end

            % Update sliders to match field values
            if app.MaxhorizontalvelSlider.Value ~= app.MaxhorizontalvelField.Value
                app.MaxhorizontalvelField.Value = app.MaxhorizontalvelSlider.Value;
            else 
            end
            if app.MaxverticalvelSlider_2.Value ~= app.MaxverticalvelField_2.Value
                app.MaxverticalvelField_2.Value = app.MaxverticalvelSlider_2.Value;
            end
            if app.JetDiameterSlider.Value ~= app.JetDiameterField.Value
                app.JetDiameterField.Value = app.JetDiameterSlider.Value;
            else
            end
            if app.zSlider.Value ~= app.zField.Value
                app.zField.Value = app.zSlider.Value;
            else 
            end

            % Update field value if check box is checked
            if app.FixedParametersCheckBox.Value == 1
                app.MaxhorizontalvelField_2.Value = app.MaxhorizontalvelField.Value;
                app.MaxhorizontalvelSlider_2.Value = app.MaxhorizontalvelSlider.Value;
                app.MaxverticalvelField.Value = app.MaxverticalvelField_2.Value;
                app.MaxverticalvelSlider.Value = app.MaxverticalvelSlider_2.Value;
                app.JetDiameterField_2.Value = app.JetDiameterField.Value;
                app.JetDiameterSlider_2.Value = app.JetDiameterSlider.Value;

                % Define parameters for vertical profile
                z = 0:0.1:1000;
                u_max = app.MaxhorizontalvelField_2.Value;
                w_max = app.MaxverticalvelField.Value;
                D = app.JetDiameterField_2.Value * 1000;
                r_max = 1.125*D;
                z_max = 0.025*D;
                r = app.rField.Value * 1000;
                ps_z = (2*r.^2./D^2).^0.35;
                gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
                de = 2.1;
                eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
                kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
                ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%                 de = 2.1;
%                 gam = 0.35;
%                 eps = 0.18;
%                 kap = 0.06;
%                 ch = 1;
                c1 = -0.16;
                c2 = 1.19;
                sf_u = u_max / (0.338 * r_max);
                sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));
    
                % Individual parts of the functions
                uz_1 = exp(-(0.7 - (2.*r.^2./D^2).^0.35).^2) + 0.18*exp(-0.06.*(r./D).^2);
                uz_2 = (z./z_max).^(c2-1) .* exp(c1.*(z./z_max).^c2);
                wz_1 = (1 + 0.7*ps_z*(0.7-ps_z)) * exp(-(0.7 - ps_z)^2) + 0.18*exp(-0.06*(r/D)^2) .* (1 - kap.*ch.*((r/D)^2));
                wz_2 = (z_max/(c1*c2)) * (exp(c1.*(z./z_max).^c2)-1);
    
                % Plot vertical profile horizontal and vertical velocities as well as magnitude
                plot(app.VerticalAxes, compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2), z, 'k', compute_uz(app, r, sf_u, uz_1, uz_2), z, 'r', compute_wz(app, sf_w, wz_1, wz_2), z, 'b');
                p = app.VerticalAxes.Children;
                p(1).LineWidth = 1.5;
                p(2).LineWidth = 1.5;
                p(3).LineWidth = 1.5;
                legend(app.VerticalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
                legend(app.VerticalAxes, 'Location', 'southoutside');

                % "Plotting options" menu check box adjustments
                if app.PlotmagnitudeCheckBox.Value == 1
                    p(3).Visible = 'on';
                else p(3).Visible = 'off';
                end
                if app.PlothorizontalvelocityCheckBox.Value == 1
                    p(2).Visible = 'on';
                else p(2).Visible = 'off';
                end
                if app.PlotverticalvelocityCheckBox.Value == 1
                    p(1).Visible = 'on';
                else p(1).Visible = 'off';
                end
            else
            end
        end

        % Value changed function: JetDiameterField_2, 
        % MaxhorizontalvelField_2, MaxverticalvelField, rField
        function FieldValue2Changed(app, event)
            % Define parameters for vertical profile
            z = 0:0.1:1000;
            u_max = app.MaxhorizontalvelField_2.Value;
            w_max = app.MaxverticalvelField.Value;
            D = app.JetDiameterField_2.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            r = app.rField.Value * 1000;
            ps_z = (2*r.^2./D^2).^0.35;
            gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
            de = 2.1;
            eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
            kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
            ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%             de = 2.1;
%             gam = 0.35;
%             eps = 0.18;
%             kap = 0.06;
%             ch = 1;
            c1 = -0.16;
            c2 = 1.19;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            uz_1 = exp(-(0.7 - (2.*r.^2./D^2).^0.35).^2) + 0.18*exp(-0.06.*(r./D).^2);
            uz_2 = (z./z_max).^(c2-1) .* exp(c1.*(z./z_max).^c2);
            wz_1 = (1 + 0.7*ps_z*(0.7-ps_z)) * exp(-(0.7 - ps_z)^2) + 0.18*exp(-0.06*(r/D)^2) .* (1 - kap.*ch.*((r/D)^2));
            wz_2 = (z_max/(c1*c2)) * (exp(c1.*(z./z_max).^c2)-1);

            % Plot vertical profile horizontal and vertical velocities as well as magnitude
            plot(app.VerticalAxes, compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2), z, 'k', compute_uz(app, r, sf_u, uz_1, uz_2), z, 'r', compute_wz(app, sf_w, wz_1, wz_2), z, 'b');
            p = app.VerticalAxes.Children;
            p(1).LineWidth = 1.5;
            p(2).LineWidth = 1.5;
            p(3).LineWidth = 1.5;
            legend(app.VerticalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
            legend(app.VerticalAxes, 'Location', 'southoutside');

            % "Plotting options" menu check box adjustments
            if app.PlotmagnitudeCheckBox.Value == 1
                p(3).Visible = 'on';
            else p(3).Visible = 'off';
            end
            if app.PlothorizontalvelocityCheckBox.Value == 1
                p(2).Visible = 'on';
            else p(2).Visible = 'off';
            end
            if app.PlotverticalvelocityCheckBox.Value == 1
                p(1).Visible = 'on';
            else p(1).Visible = 'off';
            end

            % Update field values to match sliders
            if app.MaxhorizontalvelField_2.Value ~= app.MaxhorizontalvelSlider_2.Value
                app.MaxhorizontalvelSlider_2.Value = app.MaxhorizontalvelField_2.Value;
            end
            if app.MaxverticalvelField.Value ~= app.MaxverticalvelSlider.Value
                app.MaxverticalvelSlider.Value = app.MaxverticalvelField.Value;
            else
            end
            if app.JetDiameterField_2.Value ~= app.JetDiameterSlider_2.Value
                app.JetDiameterSlider_2.Value = app.JetDiameterField_2.Value;
            else
            end
            if app.rField.Value ~= app.rSlider.Value
                app.rSlider.Value = app.rField.Value;
            else
            end
            
            % Update field value if check box is checked
            if app.FixedParametersCheckBox.Value == 1
                app.MaxhorizontalvelField.Value = app.MaxhorizontalvelField_2.Value;
                app.MaxhorizontalvelSlider.Value = app.MaxhorizontalvelSlider_2.Value;
                app.MaxverticalvelField_2.Value = app.MaxverticalvelField.Value;
                app.MaxverticalvelSlider_2.Value = app.MaxverticalvelSlider.Value;
                app.JetDiameterField.Value = app.JetDiameterField_2.Value;
                app.JetDiameterSlider.Value = app.JetDiameterSlider_2.Value;

                % Define parameters for radial profile
                r = 0:0.1:10000;
                u_max = app.MaxhorizontalvelField.Value;
                w_max = app.MaxverticalvelField_2.Value;
                D = app.JetDiameterField.Value * 1000;
                r_max = 1.125*D;
                z_max = 0.025*D;
                z = app.zField.Value;
                de = 2.1;
                gam = 0.35;
                eps = 0.18;
                kap = 0.06;
                ch = 1;
                c1 = -4.1;
                c2 = 1.1;
                ps_r = (de*r.^2./r_max^2).^gam;
                sf_u = u_max / (0.338 * r_max);
                sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));
        
                % Individual parts of the functions
                ur_1 = exp(-(2*gam - (de.*r.^2./r_max^2).^gam).^2) + eps .* exp(-kap.*(r.^2./r_max^2));
                ur_2 = (z/D)^0.1 * exp(c1*(z/D)^1.1);
                wr_1 = (1 + 2.*gam.*ps_r.*(2*gam-ps_r)) .* exp(-(2*gam - ps_r).^2) + eps .* exp(-kap.*(r.^2./r_max^2).^ch) .* (1 - kap.*ch.*(r./r_max).^(2*ch));
                wr_2 = (D/-4.51) * (exp(-4.51*(z/D).^1.1)-1);
        
                % Plot radial profile horizontal and vertical velocities as well as magnitude
                plot(app.HorizontalAxes, r, compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2), 'k', r, compute_ur(app, r, sf_u, ur_1, ur_2), 'r', r, compute_wr(app, sf_w, wr_1, wr_2), 'b');
                p = app.HorizontalAxes.Children;
                p(1).LineWidth = 1.5;
                p(2).LineWidth = 1.5;
                p(3).LineWidth = 1.5;
                legend(app.HorizontalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
                legend(app.HorizontalAxes, 'Location', 'southoutside'); 

    %             % Adjust axis for radial profile
    %             r_lim = 2 * r_max;
    %             xlim(app.HorizontalAxes, [0, (r_lim+2000)]);
    
                % "Plotting options" menu check box adjustments
                if app.PlotmagnitudeCheckBox.Value == 1
                    p(3).Visible = 'on';
                else p(3).Visible = 'off';
                end
                if app.PlothorizontalvelocityCheckBox.Value == 1
                    p(2).Visible = 'on';
                else p(2).Visible = 'off';
                end
                if app.PlotverticalvelocityCheckBox.Value == 1
                    p(1).Visible = 'on';
                else p(1).Visible = 'off';
                end    
            else
            end
        end

        % Value changed function: JetDiameterSlider_2, 
        % MaxhorizontalvelSlider_2, MaxverticalvelSlider, rSlider
        function SliderValue2Changed(app, event)
            % Define parameters for vertical profile
            z = 0:0.1:1000;
            u_max = app.MaxhorizontalvelSlider_2.Value;
            w_max = app.MaxverticalvelSlider.Value;
            D = app.JetDiameterSlider_2.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            r = app.rSlider.Value * 1000;
            ps_z = (2*r.^2./D^2).^0.35;
            gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
            de = 2.1;
            eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
            kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
            ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%             de = 2.1;
%             gam = 0.35;
%             eps = 0.18;
%             kap = 0.06;
%             ch = 1;
            c1 = -0.16;
            c2 = 1.19;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            uz_1 = exp(-(0.7 - (2.*r.^2./D^2).^0.35).^2) + 0.18*exp(-0.06.*(r./D).^2);
            uz_2 = (z./z_max).^(c2-1) .* exp(c1.*(z./z_max).^c2);
            wz_1 = (1 + 0.7*ps_z*(0.7-ps_z)) * exp(-(0.7 - ps_z)^2) + 0.18*exp(-0.06*(r/D)^2) .* (1 - kap.*ch.*((r/D)^2));
            wz_2 = (z_max/(c1*c2)) * (exp(c1.*(z./z_max).^c2)-1);

            % Plot vertical profile horizontal and vertical velocities as well as magnitude
            plot(app.VerticalAxes, compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2), z, 'k', compute_uz(app, r, sf_u, uz_1, uz_2), z, 'r', compute_wz(app, sf_w, wz_1, wz_2), z, 'b');
            p = app.VerticalAxes.Children;
            p(1).LineWidth = 1.5;
            p(2).LineWidth = 1.5;
            p(3).LineWidth = 1.5;
            legend(app.VerticalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
            legend(app.VerticalAxes, 'Location', 'southoutside');

            % "Plotting options" menu check box adjustments
            if app.PlotmagnitudeCheckBox.Value == 1
                p(3).Visible = 'on';
            else p(3).Visible = 'off';
            end
            if app.PlothorizontalvelocityCheckBox.Value == 1
                p(2).Visible = 'on';
            else p(2).Visible = 'off';
            end
            if app.PlotverticalvelocityCheckBox.Value == 1
                p(1).Visible = 'on';
            else p(1).Visible = 'off';
            end  

            % Update field values to match sliders
            if app.MaxhorizontalvelSlider_2.Value ~= app.MaxhorizontalvelField_2.Value
                app.MaxhorizontalvelField_2.Value = app.MaxhorizontalvelSlider_2.Value;
            end
            if app.MaxverticalvelSlider.Value ~= app.MaxverticalvelField.Value
                app.MaxverticalvelField.Value = app.MaxverticalvelSlider.Value;
            else
            end
            if app.JetDiameterSlider_2.Value ~= app.JetDiameterField_2.Value
                app.JetDiameterField_2.Value = app.JetDiameterSlider_2.Value;
            else
            end
            if app.rSlider.Value ~= app.rField.Value
                app.rField.Value = app.rSlider.Value;
            else
            end
            
            % Update field value if check box is checked
            if app.FixedParametersCheckBox.Value == 1
                app.MaxhorizontalvelField.Value = app.MaxhorizontalvelField_2.Value;
                app.MaxhorizontalvelSlider.Value = app.MaxhorizontalvelSlider_2.Value;
                app.MaxverticalvelField_2.Value = app.MaxverticalvelField.Value;
                app.MaxverticalvelSlider_2.Value = app.MaxverticalvelSlider.Value;
                app.JetDiameterField.Value = app.JetDiameterField_2.Value;
                app.JetDiameterSlider.Value = app.JetDiameterSlider_2.Value;

                % Define parameters for radial profile
                r = 0:0.1:10000;
                u_max = app.MaxhorizontalvelField.Value;
                w_max = app.MaxverticalvelField_2.Value;
                D = app.JetDiameterField.Value * 1000;
                r_max = 1.125*D;
                z_max = 0.025*D;
                z = app.zField.Value;
                de = 2.1;
                gam = 0.35;
                eps = 0.18;
                kap = 0.06;
                ch = 1;
                c1 = -4.1;
                c2 = 1.1;
                ps_r = (de*r.^2./r_max^2).^gam;
                sf_u = u_max / (0.338 * r_max);
                sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));
        
                % Individual parts of the functions
                ur_1 = exp(-(2*gam - (de.*r.^2./r_max^2).^gam).^2) + eps .* exp(-kap.*(r.^2./r_max^2));
                ur_2 = (z/D)^0.1 * exp(c1*(z/D)^1.1);
                wr_1 = (1 + 2.*gam.*ps_r.*(2*gam-ps_r)) .* exp(-(2*gam - ps_r).^2) + eps .* exp(-kap.*(r.^2./r_max^2).^ch) .* (1 - kap.*ch.*(r./r_max).^(2*ch));
                wr_2 = (D/-4.51) * (exp(-4.51*(z/D).^1.1)-1);
        
                % Plot radial profile horizontal and vertical velocities as well as magnitude
                plot(app.HorizontalAxes, r, compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2), 'k', r, compute_ur(app, r, sf_u, ur_1, ur_2), 'r', r, compute_wr(app, sf_w, wr_1, wr_2), 'b');
                p = app.HorizontalAxes.Children;
                p(1).LineWidth = 1.5;
                p(2).LineWidth = 1.5;
                p(3).LineWidth = 1.5;
                legend(app.HorizontalAxes, 'Magnitude', 'Horizontal Velocity', 'Vertical Velocity');
                legend(app.HorizontalAxes, 'Location', 'southoutside');             
              
    %             % Adjust axis for radial profile
    %             r_lim = 2 * r_max;
    %             xlim(app.HorizontalAxes, [0, (r_lim+2000)]);

                % "Plotting options" menu check box adjustments
                if app.PlotmagnitudeCheckBox.Value == 1
                    p(3).Visible = 'on';
                else p(3).Visible = 'off';
                end
                if app.PlothorizontalvelocityCheckBox.Value == 1
                    p(2).Visible = 'on';
                else p(2).Visible = 'off';
                end
                if app.PlotverticalvelocityCheckBox.Value == 1
                    p(1).Visible = 'on';
                else p(1).Visible = 'off';
                end    
            else
            end
        end

        % Button pushed function: ExportPlotButton
        function ExportPlotButtonPushed(app, event)
            % Create a temporary figure of the plot
            fig = figure;
            fig.Visible = 'off';
            figAxes = axes(fig);            
            allChildren = app.HorizontalAxes.XAxis.Parent.Children;
            copyobj(allChildren, figAxes)
            figAxes.XLim = app.HorizontalAxes.XLim;
            figAxes.YLim = app.HorizontalAxes.YLim;
            figAxes.ZLim = app.HorizontalAxes.ZLim;
            figAxes.DataAspectRatio = app.HorizontalAxes.DataAspectRatio;
            xlabel(figAxes, 'Distance from downdraft center [m]');
            ylabel(figAxes, 'Velocity [m/s]');
            title(figAxes, 'Radial Profile of Velocity');

            % Save plot as png file
            saveas(fig, 'OBV_AbdElaal_radial_profile_plot', 'png');

            % Delete the temporary figure
            delete(fig);
        end

        % Button pushed function: ExportDataTableButton
        function ExportDataTableButtonPushed(app, event)
            % Define parameters for radial profile
            r = 0:0.1:10000;
            u_max = app.MaxhorizontalvelField.Value;
            w_max = app.MaxverticalvelField_2.Value;
            D = app.JetDiameterField.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            z = app.zField.Value;
            de = 2.1;
            gam = 0.35;
            eps = 0.18;
            kap = 0.06;
            ch = 1;
            c1 = -4.1;
            c2 = 1.1;
            ps_r = (de*r.^2./r_max^2).^gam;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            ur_1 = exp(-(2*gam - (de.*r.^2./r_max^2).^gam).^2) + eps .* exp(-kap.*(r.^2./r_max^2));
            ur_2 = (z/D)^0.1 * exp(c1*(z/D)^1.1);
            wr_1 = (1 + 2.*gam.*ps_r.*(2*gam-ps_r)) .* exp(-(2*gam - ps_r).^2) + eps .* exp(-kap.*(r.^2./r_max^2).^ch) .* (1 - kap.*ch.*(r./r_max).^(2*ch));
            wr_2 = (D/-4.51) * (exp(-4.51*(z/D).^1.1)-1);

            % Put x axis values into an array
            r = [0.1:0.1:10000];
            % Compute velocities and insert into arrays
            m = [compute_u(app, r, sf_u, ur_1, ur_2, sf_w, wr_1, wr_2)];
            h = [compute_ur(app, r, sf_u, ur_1, ur_2)];
            v = [compute_wr(app, sf_w, wr_1, wr_2)];

            % Arrange arrays into a table
            ar = [r; m; h; v];
            t = array2table(ar.');
            t = renamevars(t, ["Var1", "Var2", "Var3", "Var4"], ["Distance from downdraft center [m]","Magnitude [m/s]","Horizontal Velocity [m/s]","Vertical Velocity [m/s]"]);
            
            % Write table into csv file
            writetable(t, 'OBV_AbdElaal_radial_profile_table.csv')
        end

        % Button pushed function: ExportPlotButton_2
        function ExportPlotButton_2Pushed(app, event)
            % Create a temporary figure of the plot
            fig = figure;
            fig.Visible = 'off';
            figAxes = axes(fig);            
            allChildren = app.VerticalAxes.XAxis.Parent.Children;
            copyobj(allChildren, figAxes)
            figAxes.XLim = app.VerticalAxes.XLim;
            figAxes.YLim = app.VerticalAxes.YLim;
            figAxes.ZLim = app.VerticalAxes.ZLim;
            figAxes.DataAspectRatio = app.VerticalAxes.DataAspectRatio;
            xlabel(figAxes, 'Velocity [m/s]');
            ylabel(figAxes, 'Height above ground [m]');
            title(figAxes, 'Vertical Profile of Velocity');

            % Save as png file
            saveas(fig, 'OBV_AbdElaal_vertical_profile_plot', 'png');

            % Delete the temporary figure
            delete(fig);
        end

        % Button pushed function: ExportDataTableButton_2
        function ExportDataTableButton_2Pushed(app, event)
            % Define parameters for vertical profile
            z = 0:0.1:1000;
            u_max = app.MaxhorizontalvelField_2.Value;
            w_max = app.MaxverticalvelField.Value;
            D = app.JetDiameterField_2.Value * 1000;
            r_max = 1.125*D;
            z_max = 0.025*D;
            r = app.rField.Value * 1000;
            ps_z = (2*r.^2./D^2).^0.35;
            gam = -227.69.*(z./D).^3+74.553.*(z./D).^2-5.0359.*(z./D)+0.4985;
            de = 2.1;
            eps = 12.938.*(z./D).^2+2.8563.*(z./D)+0.3019;
            kap = 25.466.*(z./D).^2-1.9908.*(z./D)+0.1593;
            ch = 125.38.*(z./D).^3-46.849.*(z./D).^2+0.9688.*(z./D)+0.9932;
%             de = 2.1;
%             gam = 0.35;
%             eps = 0.18;
%             kap = 0.06;
%             ch = 1;
            c1 = -0.16;
            c2 = 1.19;
            sf_u = u_max / (0.338 * r_max);
            sf_w = w_max ./ (1.919*z_max*(exp(c1.*(z./z_max).^c2)-1));

            % Individual parts of the functions
            uz_1 = exp(-(0.7 - (2.*r.^2./D^2).^0.35).^2) + 0.18*exp(-0.06.*(r./D).^2);
            uz_2 = (z./z_max).^(c2-1) .* exp(c1.*(z./z_max).^c2);
            wz_1 = (1 + 0.7*ps_z*(0.7-ps_z)) * exp(-(0.7 - ps_z)^2) + 0.18*exp(-0.06*(r/D)^2) .* (1 - kap.*ch.*((r/D)^2));
            wz_2 = (z_max/(c1*c2)) * (exp(c1.*(z./z_max).^c2)-1);

            % Put y axis values into an array
            z = [0:0.1:1000];
            % Compute velocities and insert into arrays
            m = [compute_w(app, r, sf_u, uz_1, uz_2, sf_w, wz_1, wz_2)];
            h = [compute_uz(app, r, sf_u, uz_1, uz_2)];
            v = [compute_wz(app, sf_w, wz_1, wz_2)];

            % Arrange arrays into a table
            ar = [m; h; v; z];
            t = array2table(ar.');
            t = renamevars(t, ["Var1", "Var2", "Var3", "Var4"], ["Magnitude [m/s]","Horizontal Velocity [m/s]","Vertical Velocity [m/s]","Height above ground [m]"]);

            % Write table into csv file
            writetable(t, 'OBV_AbdElaal_vertical_profile_table.csv')
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1321 571];
            app.UIFigure.Name = 'MATLAB App';

            % Create HorizontalAxes
            app.HorizontalAxes = uiaxes(app.UIFigure);
            title(app.HorizontalAxes, 'Radial Profile of Velocity')
            xlabel(app.HorizontalAxes, 'Distance from downdraft center [m]')
            ylabel(app.HorizontalAxes, 'Velocity [m/s]')
            zlabel(app.HorizontalAxes, 'Z')
            app.HorizontalAxes.Position = [343 78 289 315];

            % Create VerticalAxes
            app.VerticalAxes = uiaxes(app.UIFigure);
            title(app.VerticalAxes, 'Vertical Profile of Velocity')
            xlabel(app.VerticalAxes, 'Velocity [m/s]')
            ylabel(app.VerticalAxes, 'Height above ground [m]')
            zlabel(app.VerticalAxes, 'Z')
            app.VerticalAxes.Position = [1015 77 289 311];

            % Create AbdElaalMillsMa2013Label
            app.AbdElaalMillsMa2013Label = uilabel(app.UIFigure);
            app.AbdElaalMillsMa2013Label.FontSize = 18;
            app.AbdElaalMillsMa2013Label.FontWeight = 'bold';
            app.AbdElaalMillsMa2013Label.Position = [12 544 588 28];
            app.AbdElaalMillsMa2013Label.Text = 'Abd-Elaal, Mills, & Ma (2013)';

            % Create RadialProfileModelParametersPanel
            app.RadialProfileModelParametersPanel = uipanel(app.UIFigure);
            app.RadialProfileModelParametersPanel.Title = 'Radial Profile: Model Parameters';
            app.RadialProfileModelParametersPanel.FontWeight = 'bold';
            app.RadialProfileModelParametersPanel.Position = [20 115 307 278];

            % Create MaximumhorizontalvelocitymsLabel
            app.MaximumhorizontalvelocitymsLabel = uilabel(app.RadialProfileModelParametersPanel);
            app.MaximumhorizontalvelocitymsLabel.HorizontalAlignment = 'right';
            app.MaximumhorizontalvelocitymsLabel.Position = [21 223 188 22];
            app.MaximumhorizontalvelocitymsLabel.Text = 'Maximum horizontal velocity [m/s]';

            % Create MaxhorizontalvelSlider
            app.MaxhorizontalvelSlider = uislider(app.RadialProfileModelParametersPanel);
            app.MaxhorizontalvelSlider.Limits = [0 60];
            app.MaxhorizontalvelSlider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.MaxhorizontalvelSlider.Position = [109 211 150 3];
            app.MaxhorizontalvelSlider.Value = 40;

            % Create MaxhorizontalvelField
            app.MaxhorizontalvelField = uieditfield(app.RadialProfileModelParametersPanel, 'numeric');
            app.MaxhorizontalvelField.Limits = [0 60];
            app.MaxhorizontalvelField.ValueChangedFcn = createCallbackFcn(app, @FieldValueChanged, true);
            app.MaxhorizontalvelField.Position = [13 193 77 22];
            app.MaxhorizontalvelField.Value = 40;

            % Create JetDiameterkmLabel
            app.JetDiameterkmLabel = uilabel(app.RadialProfileModelParametersPanel);
            app.JetDiameterkmLabel.HorizontalAlignment = 'right';
            app.JetDiameterkmLabel.Position = [24 59 100 22];
            app.JetDiameterkmLabel.Text = 'Jet Diameter [km]';

            % Create JetDiameterSlider
            app.JetDiameterSlider = uislider(app.RadialProfileModelParametersPanel);
            app.JetDiameterSlider.Limits = [0 3];
            app.JetDiameterSlider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.JetDiameterSlider.Position = [112 48 150 3];
            app.JetDiameterSlider.Value = 1;

            % Create JetDiameterField
            app.JetDiameterField = uieditfield(app.RadialProfileModelParametersPanel, 'numeric');
            app.JetDiameterField.Limits = [0 3];
            app.JetDiameterField.ValueChangedFcn = createCallbackFcn(app, @FieldValueChanged, true);
            app.JetDiameterField.Position = [16 30 77 22];
            app.JetDiameterField.Value = 1;

            % Create MaximumverticalvelocitymsLabel_2
            app.MaximumverticalvelocitymsLabel_2 = uilabel(app.RadialProfileModelParametersPanel);
            app.MaximumverticalvelocitymsLabel_2.HorizontalAlignment = 'right';
            app.MaximumverticalvelocitymsLabel_2.Position = [19 143 174 22];
            app.MaximumverticalvelocitymsLabel_2.Text = 'Maximum vertical velocity [m/s]';

            % Create MaxverticalvelSlider_2
            app.MaxverticalvelSlider_2 = uislider(app.RadialProfileModelParametersPanel);
            app.MaxverticalvelSlider_2.Limits = [0 60];
            app.MaxverticalvelSlider_2.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.MaxverticalvelSlider_2.Position = [104 130 150 3];
            app.MaxverticalvelSlider_2.Value = 30;

            % Create MaxverticalvelField_2
            app.MaxverticalvelField_2 = uieditfield(app.RadialProfileModelParametersPanel, 'numeric');
            app.MaxverticalvelField_2.Limits = [0 60];
            app.MaxverticalvelField_2.ValueChangedFcn = createCallbackFcn(app, @FieldValueChanged, true);
            app.MaxverticalvelField_2.Position = [16 112 77 22];
            app.MaxverticalvelField_2.Value = 30;

            % Create ExportPlotButton
            app.ExportPlotButton = uibutton(app.UIFigure, 'push');
            app.ExportPlotButton.ButtonPushedFcn = createCallbackFcn(app, @ExportPlotButtonPushed, true);
            app.ExportPlotButton.FontWeight = 'bold';
            app.ExportPlotButton.Position = [344 27 123 31];
            app.ExportPlotButton.Text = {'Export Plot'; ''};

            % Create ExportDataTableButton
            app.ExportDataTableButton = uibutton(app.UIFigure, 'push');
            app.ExportDataTableButton.ButtonPushedFcn = createCallbackFcn(app, @ExportDataTableButtonPushed, true);
            app.ExportDataTableButton.FontWeight = 'bold';
            app.ExportDataTableButton.Position = [510 27 123 31];
            app.ExportDataTableButton.Text = {'Export Data Table'; ''};

            % Create VerticalProfileModelParametersPanel
            app.VerticalProfileModelParametersPanel = uipanel(app.UIFigure);
            app.VerticalProfileModelParametersPanel.Title = 'Vertical Profile: Model Parameters';
            app.VerticalProfileModelParametersPanel.FontWeight = 'bold';
            app.VerticalProfileModelParametersPanel.Position = [712 116 289 272];

            % Create MaximumverticalvelocitymsLabel
            app.MaximumverticalvelocitymsLabel = uilabel(app.VerticalProfileModelParametersPanel);
            app.MaximumverticalvelocitymsLabel.HorizontalAlignment = 'right';
            app.MaximumverticalvelocitymsLabel.Position = [16 141 174 22];
            app.MaximumverticalvelocitymsLabel.Text = 'Maximum vertical velocity [m/s]';

            % Create MaxverticalvelSlider
            app.MaxverticalvelSlider = uislider(app.VerticalProfileModelParametersPanel);
            app.MaxverticalvelSlider.Limits = [0 60];
            app.MaxverticalvelSlider.ValueChangedFcn = createCallbackFcn(app, @SliderValue2Changed, true);
            app.MaxverticalvelSlider.Position = [101 128 150 3];
            app.MaxverticalvelSlider.Value = 30;

            % Create MaxverticalvelField
            app.MaxverticalvelField = uieditfield(app.VerticalProfileModelParametersPanel, 'numeric');
            app.MaxverticalvelField.Limits = [0 60];
            app.MaxverticalvelField.ValueChangedFcn = createCallbackFcn(app, @FieldValue2Changed, true);
            app.MaxverticalvelField.Position = [13 110 77 22];
            app.MaxverticalvelField.Value = 30;

            % Create JetDiameterkmLabel_2
            app.JetDiameterkmLabel_2 = uilabel(app.VerticalProfileModelParametersPanel);
            app.JetDiameterkmLabel_2.HorizontalAlignment = 'right';
            app.JetDiameterkmLabel_2.Position = [21 58 100 22];
            app.JetDiameterkmLabel_2.Text = 'Jet Diameter [km]';

            % Create JetDiameterSlider_2
            app.JetDiameterSlider_2 = uislider(app.VerticalProfileModelParametersPanel);
            app.JetDiameterSlider_2.Limits = [0 3];
            app.JetDiameterSlider_2.ValueChangedFcn = createCallbackFcn(app, @SliderValue2Changed, true);
            app.JetDiameterSlider_2.Position = [109 47 150 3];
            app.JetDiameterSlider_2.Value = 1;

            % Create JetDiameterField_2
            app.JetDiameterField_2 = uieditfield(app.VerticalProfileModelParametersPanel, 'numeric');
            app.JetDiameterField_2.Limits = [0 3];
            app.JetDiameterField_2.ValueChangedFcn = createCallbackFcn(app, @FieldValue2Changed, true);
            app.JetDiameterField_2.Position = [13 29 77 22];
            app.JetDiameterField_2.Value = 1;

            % Create MaximumhorizontalvelocitymsLabel_2
            app.MaximumhorizontalvelocitymsLabel_2 = uilabel(app.VerticalProfileModelParametersPanel);
            app.MaximumhorizontalvelocitymsLabel_2.HorizontalAlignment = 'right';
            app.MaximumhorizontalvelocitymsLabel_2.Position = [21 222 188 22];
            app.MaximumhorizontalvelocitymsLabel_2.Text = 'Maximum horizontal velocity [m/s]';

            % Create MaxhorizontalvelSlider_2
            app.MaxhorizontalvelSlider_2 = uislider(app.VerticalProfileModelParametersPanel);
            app.MaxhorizontalvelSlider_2.Limits = [0 60];
            app.MaxhorizontalvelSlider_2.ValueChangedFcn = createCallbackFcn(app, @SliderValue2Changed, true);
            app.MaxhorizontalvelSlider_2.Position = [109 210 150 3];
            app.MaxhorizontalvelSlider_2.Value = 40;

            % Create MaxhorizontalvelField_2
            app.MaxhorizontalvelField_2 = uieditfield(app.VerticalProfileModelParametersPanel, 'numeric');
            app.MaxhorizontalvelField_2.Limits = [0 60];
            app.MaxhorizontalvelField_2.ValueChangedFcn = createCallbackFcn(app, @FieldValue2Changed, true);
            app.MaxhorizontalvelField_2.Position = [13 192 77 22];
            app.MaxhorizontalvelField_2.Value = 40;

            % Create ExportPlotButton_2
            app.ExportPlotButton_2 = uibutton(app.UIFigure, 'push');
            app.ExportPlotButton_2.ButtonPushedFcn = createCallbackFcn(app, @ExportPlotButton_2Pushed, true);
            app.ExportPlotButton_2.FontWeight = 'bold';
            app.ExportPlotButton_2.Position = [1015 27 123 31];
            app.ExportPlotButton_2.Text = {'Export Plot'; ''};

            % Create ExportDataTableButton_2
            app.ExportDataTableButton_2 = uibutton(app.UIFigure, 'push');
            app.ExportDataTableButton_2.ButtonPushedFcn = createCallbackFcn(app, @ExportDataTableButton_2Pushed, true);
            app.ExportDataTableButton_2.FontWeight = 'bold';
            app.ExportDataTableButton_2.Position = [1180 27 123 31];
            app.ExportDataTableButton_2.Text = {'Export Data Table'; ''};

            % Create Panel_2
            app.Panel_2 = uipanel(app.UIFigure);
            app.Panel_2.Position = [712 34 289 83];

            % Create zmLabel_2
            app.zmLabel_2 = uilabel(app.Panel_2);
            app.zmLabel_2.HorizontalAlignment = 'right';
            app.zmLabel_2.Position = [11 53 202 22];
            app.zmLabel_2.Text = 'Distance from downdraft center [km]';

            % Create rField
            app.rField = uieditfield(app.Panel_2, 'numeric');
            app.rField.Limits = [0 3];
            app.rField.ValueChangedFcn = createCallbackFcn(app, @FieldValue2Changed, true);
            app.rField.Position = [13 24 77 22];
            app.rField.Value = 1;

            % Create rSlider
            app.rSlider = uislider(app.Panel_2);
            app.rSlider.Limits = [0 3];
            app.rSlider.ValueChangedFcn = createCallbackFcn(app, @SliderValue2Changed, true);
            app.rSlider.Position = [106 43 150 3];
            app.rSlider.Value = 1;

            % Create Panel
            app.Panel = uipanel(app.UIFigure);
            app.Panel.Position = [20 34 307 83];

            % Create zmLabel
            app.zmLabel = uilabel(app.Panel);
            app.zmLabel.HorizontalAlignment = 'right';
            app.zmLabel.Position = [23 44 138 38];
            app.zmLabel.Text = 'Height above ground [m]';

            % Create zSlider
            app.zSlider = uislider(app.Panel);
            app.zSlider.Limits = [0 1000];
            app.zSlider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.zSlider.Position = [103 43 150 3];
            app.zSlider.Value = 100;

            % Create zField
            app.zField = uieditfield(app.Panel, 'numeric');
            app.zField.Limits = [0 1000];
            app.zField.ValueChangedFcn = createCallbackFcn(app, @FieldValueChanged, true);
            app.zField.Position = [15 24 77 22];
            app.zField.Value = 100;

            % Create PlottingOptionsPanel
            app.PlottingOptionsPanel = uipanel(app.UIFigure);
            app.PlottingOptionsPanel.Title = 'Plotting Options';
            app.PlottingOptionsPanel.FontWeight = 'bold';
            app.PlottingOptionsPanel.Position = [20 463 504 61];

            % Create PlotmagnitudeCheckBox
            app.PlotmagnitudeCheckBox = uicheckbox(app.PlottingOptionsPanel);
            app.PlotmagnitudeCheckBox.Text = 'Plot magnitude';
            app.PlotmagnitudeCheckBox.Position = [13 10 103 22];

            % Create PlothorizontalvelocityCheckBox
            app.PlothorizontalvelocityCheckBox = uicheckbox(app.PlottingOptionsPanel);
            app.PlothorizontalvelocityCheckBox.Text = 'Plot horizontal velocity';
            app.PlothorizontalvelocityCheckBox.Position = [165 10 143 22];

            % Create PlotverticalvelocityCheckBox
            app.PlotverticalvelocityCheckBox = uicheckbox(app.PlottingOptionsPanel);
            app.PlotverticalvelocityCheckBox.Text = 'Plot vertical velocity';
            app.PlotverticalvelocityCheckBox.Position = [357 10 129 22];

            % Create Panel_3
            app.Panel_3 = uipanel(app.UIFigure);
            app.Panel_3.Position = [20 426 504 38];

            % Create FixedParametersCheckBox
            app.FixedParametersCheckBox = uicheckbox(app.Panel_3);
            app.FixedParametersCheckBox.Text = 'Keep model parameters fixed';
            app.FixedParametersCheckBox.Position = [13 8 180 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = OBV_AbdElaal_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end