% processExptData.m
%
% Function to take output from preprocessUserDaq.m and extract
%  appropriately scaled and named data ephys(voltage, current,
%  scaled out, gain, mode, freq) and g4 display panels (x position)
%  and fictrac (x position, y position, heading)
%
% INPUTS:
%   daqData - data collected on DAQ, with fields labeled
%   daqOutput - signal output on DAQ during experiment, with fields labeled
%   daqTime - timing vector for daqData and daqOutput
%   inputParams - input parameters from trial function (e.g. ephysRecording)
%   settings - settings struct from ephysSettings
%
% OUTPUTS:
%   exptData - struct of appropriately scaled ephys and/or behavior data
%   ephysMeta - struct of ephys metadatam from decoding telegraphed output,
%       trial parameters
%
% Created:  04/03/2021 - MC
%           11/04/2021 - MC removed g3, updated to g4 display
%           11/10/2021 - MC fixed ball width
%           11/15/2021 - MC rotated t, resampled fictrac to match DAQ
%           01/04/2022 - MC fixed resampling error w/fictrac data
%

function [exptData, exptMeta] = processExptData(daqData, daqOutput, daqTime, inputParams, settings)

%% meta data

% transfer meta data
exptMeta.exptCond = inputParams.exptCond;
exptMeta.startTimeStamp = inputParams.startTimeStamp;

% if input resistance was calculated for said trial
if isfield(inputParams,'inputResistance')
    exptMeta.inputResistance = inputParams.inputResistance;
end

% time points of recording
exptData.t = daqTime';
ln = length(exptData.t);


%% ephys data
    % check if ephys data was collected, if so, process
    if contains(inputParams.exptCond,'ephys','IgnoreCase',true)
        
        % decode telegraphed output
        exptMeta.gain = decodeTelegraphedOutput(daqData.ampGain, 'gain');
        exptMeta.freq = decodeTelegraphedOutput(daqData.ampFreq, 'freq');
        exptMeta.mode = decodeTelegraphedOutput(daqData.ampMode, 'mode');

        % process non-scaled output
        % voltage in mV
        exptData.voltage = settings.Vm.softGain .* daqData.amp10Vm;
        % current in pA
        exptData.current = settings.I.softGain .* daqData.ampI;

        % process scaled output
        switch exptMeta.mode
            case {'Track','V-Clamp'} % voltage clamp, scaled out is current
                % I = alpha * beta mV/pA (1000 is to convert from V to mV)
                exptMeta.softGain = 1000 / ...
                    (exptMeta.gain * settings.amp.beta);
                % scaled out is current, in pA
                exptData.scaledCurrent = exptMeta.softGain .* ...
                    daqData.ampScaledOut;
                % current clamp, scaled out is voltage
            case {'I=0','I-Clamp','I-Clamp Fast'}
                % V = alpha mV / mV (1000 for V to mV)
                exptMeta.softGain = 1000 / exptMeta.gain;
                % scaled out is voltage, in mV
                exptData.scaledVoltage = exptMeta.softGain .* ...
                    daqData.ampScaledOut;
        end
    end
    
%% output data
    % check if output (current inj) data was collected, if so, process
    if contains(inputParams.exptCond,'inj','IgnoreCase',true)
        
        % convert from DAQ output back to target current
        exptData.iInj = daqOutput.ampExtCmdIn/settings.VOut.IConvFactor;
    end


%% g4 panel data
    % check if panels were used, if so, process
    if contains(inputParams.exptCond,'g4','IgnoreCase',true)
        
        %load g4 settings
        userSettings 
        xpixel = NumofColumns *16;
        
        %pull xpos data
        g4xpos_px = daqData.g4display_xpos;
        %convert from pixels to degrees
        g4xpos_deg = ((double(g4xpos_px)./xpixel)* 360);
        
        exptData.g4displayXPos = g4xpos_deg;

    end
    
%% behavior data
    % check if behavior data was collected, if so, process
    if contains(inputParams.exptCond,'fictrac','IgnoreCase',true)
        
        %note:
        %X is forward
        %Y is side-to-side (roll)
        %heading is angular (yaw)
        
        fictrac_rate = 50; %set fictrac acquisition rate
        ballRadius = 9/2; %set ball radius in mm

        % pull preprocessed fictrac data in order to processes together
        positionVoltage = [daqData.ficTracHeading, daqData.ficTracIntX,...
            daqData.ficTracIntY];

        % 1)Tranform signal from voltage to radians for unwrapping
        positionRad = positionVoltage*(2*pi)./10;

        % 2)Unwrap
        positionRad_uw = unwrap(positionRad);

        % 3)Downsample the position data to match FicTrac's output
        positionRad_uw_ds = resample(positionRad_uw,(fictrac_rate/2),settings.bob.sampRate);

        % 4)Smooth the data
        positionRad_uw_ds_sm = smoothdata(positionRad_uw_ds,'rlowess',25);

        % 5)Transform to useful systems
        positionUnit(:,1) = rad2deg(positionRad_uw_ds_sm(:,1)); %degrees for yaw (0-360)
        positionUnit(:,2:3) = positionRad_uw_ds_sm(:,2:3) .* ballRadius; %mm for x/y (0-2pi*r)

        % 6)Take the derivative (must be done one at a time)
        velocity(:,1) = gradient(positionUnit(:,1)).*(fictrac_rate/2);
        velocity(:,2) = gradient(positionUnit(:,2)).*(fictrac_rate/2);
        velocity(:,3) = gradient(positionUnit(:,3)).*(fictrac_rate/2);

        % 7)OPTIONAL: Remove extreme values (must be done one at a time)
        if 0
            velocity_interpolated = zeros(length(velocity),3);
            for v=1:3
                
                % pull velocity signal one at a time
                vel_ex = velocity(:,v);
                
                % 7)Calculate the distribution and take away values that are below 2.5% and above 97.5%
                percentilelow = prctile(vel_ex,2.5);
                percentilehigh = prctile(vel_ex,97.5);
                boundedVelocity = vel_ex;
                boundedVelocity(vel_ex<percentilelow | vel_ex>percentilehigh) = NaN;
                
                % 8)Linearly interpolate to replace the NaNs with values.
                [pointsVectorV] = find(~isnan(boundedVelocity));
                valuesVectorV = boundedVelocity(pointsVectorV);
                xiV = 1:length(boundedVelocity);
                
                % load into table
                velocity_interpolated(:,v) = interp1(pointsVectorV,valuesVectorV,xiV);
            end
        else
            velocity_interpolated = velocity;
        end

        % 8)OPTIONAL: Smooth
        if 1
            velocity_sm = smoothdata(velocity_interpolated,'rlowess',15);
        else
            velocity_sm = velocity_interpolated
        end
        
        % 9)Resample to match DAQ
        % add capps to avoid end resampling error
        cap = 10;
        velocity_sm_cap = [repmat(velocity_sm(1,:),cap,1); velocity_sm; repmat(velocity_sm(end,:),cap,1)];
        position_sm_cap = [repmat(positionUnit(1,:),cap,1); positionUnit; repmat(positionUnit(end,:),cap,1)];
        % resample
        velocity_rs_cap = resample(velocity_sm_cap,settings.bob.sampRate,(50/2),3,10);
        positionUnit_rs_cap = resample(position_sm_cap,settings.bob.sampRate,(50/2),3,10);
        % remove caps
        rsFactor = cap * settings.bob.sampRate/(50/2);
        velocity_rs = velocity_rs_cap;
        velocity_rs(1:rsFactor,:) = [];
        velocity_rs(end-rsFactor+1:end,:) = [];
        positionUnit_rs = positionUnit_rs_cap;
        positionUnit_rs(1:rsFactor,:) = [];
        positionUnit_rs(end-rsFactor+1:end,:) = [];


        %Assign output structure
        exptData.headingPosition = positionUnit_rs(1:ln,1); %heading position in degrees
        exptData.angularVelocity = velocity_rs(1:ln,1); %angular velcoity in degrees/s

        exptData.XPosition = positionUnit_rs(1:ln,2); %x position in mm
        exptData.forwardVelocity = velocity_rs(1:ln,2); %forward velocity in mm/s

        exptData.YPosition = positionUnit_rs(1:ln,3); %y position in mm
        exptData.sidewaysVelocity = velocity_rs(1:ln,3); %sideways velocity in mm/s
        
        %optional plot for checking analysis
        if 0
            for idx = 2
            figure
            
            subplot(9,1,1),plot(positionVoltage(:,idx)),title('voltage'), xlim([0 length(positionVoltage)])
            subplot(9,1,2),plot(positionRad(:,idx)),title('radians'), xlim([0 length(positionRad)])
            subplot(9,1,3),plot(positionRad_uw(:,idx)),title('radians unwrapped'), xlim([0 length(positionRad_uw)])
            subplot(9,1,4),plot(positionRad_uw_ds_sm(:,idx)),title('radians downsampled/smoothed'), xlim([0 length(positionRad_uw_ds_sm)])
            subplot(9,1,5),plot(positionUnit(:,idx)),title('units of interest'), xlim([0 length(positionUnit)])
            subplot(9,1,6),plot(velocity(:,idx)),title('velocity'), xlim([0 length(velocity)])
            subplot(9,1,7),plot(velocity_interpolated(:,idx)),title('velocity int'), xlim([0 length(velocity_interpolated)])
            subplot(9,1,8),plot(velocity_sm(:,idx)),title('velocity sm'), xlim([0 length(velocity_sm)])
            subplot(9,1,9),plot(velocity_rs(:,idx)),title('velocity resampled'), xlim([0 length(velocity_rs)])
            
            end
        end
       
    end

    

end