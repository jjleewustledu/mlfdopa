classdef Ki_Criswell  
    %% KI_CRISWELL estimates Ki by Patlak's method from TACs from various regions of interest, one of which is 
    %  specified to be a reference region.  
    %  Created 2018 Dec 11.  Developed on Matlab 9.4.0.813654 (R2018a).  Copyright 2018 John Joowon Lee.

    properties
        doplots
        suffix
        writecsv
    end
    
    properties (Dependent)
        csv
        estimates
        Ki
        CIupper
        CIlower
        pValue
        Rsquared
        Nregions
        sampleStartTime % min
        sampleEndTime   % min
    end
    
    methods (Static)
        function [T,this] = Kocc(varargin)
            %% KOCC estimates Ki by Patlak's method from TACs from various regions of interest with respect to the TAC
            %  from the occipital region.  KI_CRISWELL expects the TACs to be organized in a csv file as follows.
            %  - First row contains labels:
            %    patid,frame,Length(sec),Midpoint(sec),occip, ...
            %    L_caudate,R_caudate,L_ant_putamen,R_ant_putamen,L_pos_putamen,R_pos_putamen,L_thalamus,R_thalamus,...
            %    L_pallidum,R_pallidum,L_VDC,R_VDC.
            %  - Subsequent rows contain comma separated numerical values.
            %
            %  The implementation generates Patlak quantities including the ratio of tracer activities 
            %  compared to reference activities and cumulative time integrals of the reference activities;
            %  it then fits a linear regression with 5-fold cross-validation.
            %  See also:  Phelps, M.  PET:  Molecular Imaging and Biological Applications.  (2004)  pp. 165--169 and refs.
            %
            %  @param (required) csv is the name of the csv file.
            %  @param named sampleStartTime is the start time of activity samples for a Patlak model (min).
            %         The requested start time will not precede the peak of the reference activity curve.
            %  @param named sampleEndTime is the end time of activity samples for a Patlak model (min).
            %         The requested end time will not exceed the end time for the data.
            %  @param named writecsv is logical; default == true; if true write results in a new csv named
            %         using the fileprefix of param csv + param suffix + '.csv'.
            %  @param named suffix is char; default == '_Kocc'.
            %  @param named refreg is a string for the reference which matches one of the regions listed in the first 
            %         row of labels; default is 'occip'.
            %  @param named doplots is logical; default == false;
            %  @return Kocc and statistics estimated by linear regression according to Patlak's method (1/min).
            %  @return this object which may be queried for diagnostic information.
            
            if (isempty(varargin))
                if (8 == exist('Ki_Criswell', 'class'))
                    disp(evalc('help Ki_Criswell.Kocc'))
                elseif (8 == exist('mlfdopa.Ki_Criswell', 'class'))
                    disp(evalc('help mlfdopa.Ki_Criswell.Kocc'))
                end    
                error('mlfdopa:RuntimeWarning', ...
                    'Ki_Criswell.Kocc:  needs more information regarding what you wish to do.');               
            end
            if (~ischar(varargin{1}))
                fprintf('Ki_Criswell.Kocc:  doesn''t support specifying CSV files expressed as type %s\n', class(varargin{1}));
                error('mlfdopa:Value', ...
                    'Ki_Criswell.Kocc:  needs the name of an existing CSV files expressed as char.');
            end
            if (2 ~= exist(varargin{1}, 'file'))
                fprintf('Ki_Criswell.Kocc:  could not find a file named %s.  ', varargin{1}); 
                fprintf('Your current working directory is %s.\n', pwd);
                error('mlfdopa:FileNotFoundError', ...
                    'Please check the directory and name of the file you need.');
            end
            if (8 == exist('Ki_Criswell', 'class'))
                this = Ki_Criswell(varargin{:}); 
            elseif (8 == exist('mlfdopa.Ki_Criswell', 'class'))
                this = mlfdopa.Ki_Criswell(varargin{:}); 
            else
                disp('Please see more about <a href = "https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html">The Matlab Search Path</a>');
                error('mlfdopa:ValueError', ...
                    'Ki_Criswell.Kocc:  the requested class, Ki_Criswell, is not in Matlab''s path');
            end 
            try   
                T = this.writetable;
                if (this.doplots)
                    this.plotTacs;
                    this.plotIntermed;
                    this.plotPatlak;     
                end
            catch ME  
                dispexcept(ME, 'mlfdopa:RuntimeWarning', 'Ki_Criswell.Kocc');
            end
        end
    end
    
    methods 
        function this = Ki_Criswell(varargin)
            %% KI_CRISWELL estimates Ki by Patlak's method from TACs from various regions of interest, one of which is 
            %  specified to be a reference region.  KI_CRISWELL expects the TACs to be organized in a csv file as follows.
            %  - First row contains labels:
            %    patid,frame,Length(sec),Midpoint(sec),occip, ...
            %    L_caudate,R_caudate,L_ant_putamen,R_ant_putamen,L_pos_putamen,R_pos_putamen,L_thalamus,R_thalamus,...
            %    L_pallidum,R_pallidum,L_VDC,R_VDC.
            %  - Subsequent rows contain comma separated numerical values.
            %
            %  The implementation generates Patlak quantities including the ratio of tracer activities 
            %  compared to reference activities and cumulative time integrals of the reference activities;
            %  it then fits a linear regression with 5-fold cross-validation.
            %  See also:  Phelps, M.  PET:  Molecular Imaging and Biological Applications.  (2004)  pp. 165--169 and refs.
            %
            %  @param (required) csv is the name of the csv file.
            %  @param named sampleStartTime is the start time of activity samples for a Patlak model (min).
            %         The requested start time will not precede the peak of the reference activity curve.
            %  @param named sampleEndTime is the end time of activity samples for a Patlak model (min).
            %         The requested end time will not exceed the end time for the data.
            %  @param named writecsv is logical; default == true; if true write results in a new csv named
            %         using the fileprefix of param csv + param suffix + '.csv'.
            %  @param named suffix is char; default == '_Kocc'.
            %  @param named refreg is a string for the reference which matches one of the regions listed in the first 
            %         row of labels; default is 'occip'.
            %  @param named doplots is logical; default == false;
            %  @return this object which may be queried for diagnostic information.
            
            ip = inputParser;
            addRequired(ip, 'csv', @(x) 2 == exist(x, 'file') || 2 == exist([x '.csv'], 'file'));
            addParameter(ip, 'sampleStartTime', 0, @isnumeric);
            addParameter(ip, 'sampleEndTime', Inf, @isnumeric);
            addParameter(ip, 'writecsv', true, @islogical);
            addParameter(ip, 'suffix', '_Kocc', @ischar);
            addParameter(ip, 'refreg', 'occip', @ischar);
            addParameter(ip, 'doplots', false, @islogical);
            parse(ip, varargin{:});
            this.csv_ = ip.Results.csv;
            this.writecsv = ip.Results.writecsv;
            this.suffix  = ip.Results.suffix;
            this.doplots = ip.Results.doplots;
            
            % configure table objects
            this.tbl     = readtable(ip.Results.csv);
            this.t       = this.tbl.Midpoint_sec_;
            this.refreg  = strrep(ip.Results.refreg, '(', '_');
            this.refreg  = strrep(this.refreg,            ')', '_');
            this.regions = this.setdiff( ...
                this.tbl.Properties.VariableNames, [{'patid' 'frame' 'Length_sec_' 'Midpoint_sec_'} this.refreg]);    
            
            % estimate paramters
            this.frame0_ = max(this.findFrameOfMaxCp, this.findFrameOf(ip.Results.sampleStartTime*60));
            this.framef_ = min(this.findFrameOf(ip.Results.sampleEndTime*60), length(this.t));            
            this.patlakt = this.int_Cp_over_Cp;
            ratio = this.TAC_over_Cp;
            sz = size(ratio);
            for k = 1:sz(2)
                e(k) = this.estimateKi(ratio(:,k)); %#ok<AGROW>
            end
            this.estimates_ = e;
        end
        
        %% GET
        
        function g = get.csv(this)
            g = this.csv_;
        end
        function g = get.estimates(this)
            g = this.estimates_;
        end
        function g = get.Ki(this)
            g = struct2table(this.estimates);
            g = g.Ki;
        end
        function g = get.CIupper(this)
            g = struct2table(this.estimates);
            g = g.CIupper;
        end
        function g = get.CIlower(this)
            g = struct2table(this.estimates);
            g = g.CIlower;
        end
        function g = get.pValue(this)
            g = struct2table(this.estimates);
            g = g.pValue;
        end
        function g = get.Rsquared(this)
            g = struct2table(this.estimates);
            g = g.Rsquared;
        end
        function g = get.Nregions(this)
            g = size(this.regions, 2);
        end 
        function g = get.sampleStartTime(this)
            g = floor(this.t(this.frame0_)/60);
        end 
        function g = get.sampleEndTime(this)
            g = floor(this.t(this.framef_)/60);
        end 
        
        %%
        
        function plotIntermed(this)
            figure
            plot(this.t/60, this.TAC_over_Cp);
            title(['Ki_Criswell.plotIntermed:  ' ...
                  sprintf('sampling %g-%g min', this.sampleStartTime, this.sampleEndTime)], 'Interpreter', 'none');
            xlabel('time / min', 'FontSize', 14);
            ylabel(sprintf('A_{TAC}(t) / A_{%s}(t)', this.refreg), 'FontSize', 14)
            figure
            title(['Ki_Criswell.plotIntermed:  ' ...
                  sprintf('sampling %g-%g min', this.sampleStartTime, this.sampleEndTime)], 'Interpreter', 'none');
            plot(this.t/60, this.int_Cp_over_Cp, '-+');
            xlabel('time / min', 'FontSize', 14);
            ylabel(sprintf('\\int^t dt'' A_{%s}(t'') / A_{%s}(t)', this.refreg, this.refreg), 'FontSize', 14)            
        end
        function plotPatlak(this)
            figure
            plot(this.int_Cp_over_Cp/60, this.TAC_over_Cp);
            idx0 = this.findFrameOf(this.sampleStartTime*60);
            idxf = this.findFrameOf(this.sampleEndTime*60);
            title(['Ki_Criswell.plotPatlak:  ' ...
                  sprintf('sampling cumulative %g-%g min', floor(this.patlakt(idx0)/60), floor(this.patlakt(idxf)/60))], 'Interpreter', 'none');
            xlabel(sprintf('\\int^t dt'' A_{%s}(t'') / A_{%s}(t) / min', this.refreg, this.refreg), 'FontSize', 14)   
            ylabel('A_{TAC}(t) / A_{ref}(t)', 'FontSize', 14)
            legend(this.regions, 'Interpreter', 'none', 'Location', 'east');
        end
        function plotTacs(this)
            figure
            plot(this.t/60, table2array(this.tbl(:, this.refreg)), '-+', 'LineWidth', 3);
            title(['Ki_Criswell.plotTacs:  ' ...
                  sprintf('sampling %g-%g min', this.sampleStartTime, this.sampleEndTime)], 'Interpreter', 'none');
            hold on
            plot(this.t/60, table2array(this.tbl(:, this.regions)));
            xlabel('time / min', 'FontSize', 14);
            ylabel('specific activity / (Bq/mL)', 'FontSize', 14)
            legend([this.refreg this.regions], 'Interpreter', 'none', 'Location', 'east');
            hold off
        end
        function T = writetable(this)
            c = squeeze(struct2cell(this.estimates));
            T = table(cell2mat(c(:,1)), ...
                    'VariableNames', {this.regions{1}}, ...
                    'RowNames', {'Kocc' 'upper confidence' 'lower confidence' 'p-value' 'R-squared'});
            for icol = 2:size(c,2)
                T = horzcat(T, table(cell2mat(c(:,icol)), 'VariableNames', {this.regions{icol}})); %#ok<AGROW>
            end
            [pth,fp] = fileparts(this.csv);
            if (this.writecsv)
                writetable(T, fullfile(pth, [fp this.suffix '.csv']));
            end
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        csv_
        frame0_
        framef_
        estimates_
        patlakt
        refreg
        regions
        tbl
        t
    end
    
    methods (Access = protected)
        function e = estimateKi(this, r_)
            %  @param r_ is the Patlak activities ratio.
            
            t_ = this.patlakt;
            f0 = this.frame0_; 
            ff = this.framef_; 
            training = table(t_(f0:ff), r_(f0:ff), 'VariableNames', {'t', 'ratio'});
            model = this.trainRegressionModel(training);
            
            coefCI_    = model.LinearModel.coefCI * 60; % 1/min
            e.Ki       = table2array(model.LinearModel.Coefficients('t', 'Estimate')) * 60; % 1/min
            e.CIupper  = coefCI_(2,2);
            e.CIlower  = coefCI_(2,1);
            e.pValue   = table2array(model.LinearModel.Coefficients('t', 'pValue'));
            e.Rsquared = model.LinearModel.Rsquared.Ordinary;
        end
        function fr = findFrameOf(this, timepoint)
            %% of wall-clock time
            
            if (timepoint >= this.t(end))
                fr = length(this.t);
                return
            end
            [~,fr] = max(this.t > timepoint);
        end
        function fr = findFrameOfMaxCp(this)
            %% argmax_i C_{\text{plasma}}(t_i)

            cp     = table2array(this.tbl(:, this.refreg));  
            [~,fr] = max(cp);
        end
        function patlakt = int_Cp_over_Cp(this)
            %% Patlak time:  $\int^t dt' C_{\text{plasma}}(t') / C_{\text{plasma}}(t)$
            %  @return ptimes is double array of size N_times (x) 1.

            cp      = table2array(this.tbl(:, this.refreg));        
            patlakt = cumtrapz(this.t, cp) ./ cp;
        end
        function regs1 = setdiff(~, regs0, toexcl)
            regs1 = {};
            for r = 1:length(regs0)
                if (~contains(regs0{r}, toexcl))
                    regs1 = [regs1 regs0{r}]; %#ok<AGROW>
                end
            end            
        end
        function ratio = TAC_over_Cp(this)
            %% Ratio of dynamic activities:  $C_{\text{TAC}}(t) / C_{\text{plasma}}(t)$.  
            %  @return ratio is double array of size N_times (x) N_regions.

            tac   = table2array(this.tbl(:, this.regions));
            cp    = table2array(this.tbl(:, this.refreg));
            ratio = zeros(size(cp));
            for r = 1:size(tac, 2)
                ratio(:,r) = tac(:,r) ./ cp;
            end
        end
        
        %% MATLAB generated

        function [trainedModel, validationRMSE] = trainRegressionModel(~, trainingData)
            % [trainedModel, validationRMSE] = trainRegressionModel(trainingData)
            % returns a trained regression model and its RMSE. This code recreates the
            % model trained in Regression Learner app. Use the generated code to
            % automate training the same model with new data, or to learn how to
            % programmatically train models.
            %
            %  Input:
            %      trainingData: a table containing the same predictor and response
            %       columns as imported into the app.
            %
            %  Output:
            %      trainedModel: a struct containing the trained regression model. The
            %       struct contains various fields with information about the trained
            %       model.
            %
            %      trainedModel.predictFcn: a function to make predictions on new data.
            %
            %      validationRMSE: a double containing the RMSE. In the app, the
            %       History list displays the RMSE for each model.
            %
            % Use the code to train the model with new data. To retrain your model,
            % call the function from the command line with your original data or new
            % data as the input argument trainingData.
            %
            % For example, to retrain a regression model trained with the original data
            % set T, enter:
            %   [trainedModel, validationRMSE] = trainRegressionModel(T)
            %
            % To make predictions with the returned 'trainedModel' on new data T2, use
            %   yfit = trainedModel.predictFcn(T2)
            %
            % T2 must be a table containing at least the same predictor columns as used
            % during training. For details, enter:
            %   trainedModel.HowToPredict

            % Auto-generated by MATLAB on 12-Dec-2018 20:08:52


            % Extract predictors and response
            % This code processes the data into the right shape for training the
            % model.
            inputTable = trainingData;
            predictorNames = {'t'};
            predictors = inputTable(:, predictorNames);
            response = inputTable.ratio;
            isCategoricalPredictor = [false]; %#ok<*NBRAK,*NASGU>

            % Train a regression model
            % This code specifies all the model options and trains the model.
            concatenatedPredictorsAndResponse = predictors;
            concatenatedPredictorsAndResponse.ratio = response;
            linearModel = fitlm(...
                concatenatedPredictorsAndResponse, ...
                'linear', ...
                'RobustOpts', 'off');

            % Create the result struct with predict function
            predictorExtractionFcn = @(t) t(:, predictorNames);
            linearModelPredictFcn = @(x) predict(linearModel, x);
            trainedModel.predictFcn = @(x) linearModelPredictFcn(predictorExtractionFcn(x));

            % Add additional fields to the result struct
            trainedModel.RequiredVariables = {'t'};
            trainedModel.LinearModel = linearModel;
            trainedModel.About = 'This struct is a trained model exported from Regression Learner R2018a.';
            trainedModel.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

            % Extract predictors and response
            % This code processes the data into the right shape for training the
            % model.
            inputTable = trainingData;
            predictorNames = {'t'};
            predictors = inputTable(:, predictorNames);
            response = inputTable.ratio;
            isCategoricalPredictor = [false];

            % Perform cross-validation
            KFolds = 5;
            cvp = cvpartition(size(response, 1), 'KFold', KFolds);
            % Initialize the predictions to the proper sizes
            validationPredictions = response;
            for fold = 1:KFolds
                trainingPredictors = predictors(cvp.training(fold), :);
                trainingResponse = response(cvp.training(fold), :);
                foldIsCategoricalPredictor = isCategoricalPredictor;

                % Train a regression model
                % This code specifies all the model options and trains the model.
                concatenatedPredictorsAndResponse = trainingPredictors;
                concatenatedPredictorsAndResponse.ratio = trainingResponse;
                linearModel = fitlm(...
                    concatenatedPredictorsAndResponse, ...
                    'linear', ...
                    'RobustOpts', 'off');

                % Create the result struct with predict function
                linearModelPredictFcn = @(x) predict(linearModel, x);
                validationPredictFcn = @(x) linearModelPredictFcn(x);

                % Add additional fields to the result struct

                % Compute validation predictions
                validationPredictors = predictors(cvp.test(fold), :);
                foldPredictions = validationPredictFcn(validationPredictors);

                % Store predictions in the original order
                validationPredictions(cvp.test(fold), :) = foldPredictions;
            end

            % Compute validation RMSE
            isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
            validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));

        end
    end
end