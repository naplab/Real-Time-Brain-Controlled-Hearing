%% About

% This script is used to load a trained model for a subject and run
% real-time inference with it.

% The model reconstructs the envelope of the attended speech based on the
% subject's invasive brain signals and compares it with the envelopes of
% the conversations in the scene. The conversation that yields the highest
% match is decoded as the attended conversation.

%% Prepare the workspace

clc;
close all;
clear all;

%% Change to current directory 

cd C:\Users\Documents\'NimaMesgarani Lab';

% tmp = matlab.desktop.edit.getActive;
% ccd(fileparts(tmp.Filename));

%% Add path to TDT's MATLAB functions + custom functions

addpath(genpath("C:\TDT\TDTMatlabSDK"));
addpath(genpath(".\Functions"));

%% Create UDP Port

% warning("Testing mode!");
% u = udpport;

u = udpport("byte", "IPV4", "LocalHost", '192.168.1.2', "LocalPort", 1112);

%% Run the experiment file with commands sent through SynapseAPI

syn = SynapseAPI('localhost');

%% Setup the APIStreamer
% Make sure you have the online file opened in Synapse!

syn.setModeStr('Record');
s = APIStreamer('GIZMO', 'APIStreamerMC', 'HISTORY', 1, 'CALLBACK', @my_api_callback);

%% Design filters for HG Env

fs3 = 400; % Hz
notchFreq = 120; % Hz
[b,a] = fir1(100,2*[70 150]/fs3);
[bl, al] = fir1(100, 2*20/fs3); % Changed to 20 Hz
[bN, aN] = butter(5, [notchFreq - 1 notchFreq + 1]/(fs3/2), "stop");

%% Design filters for LF Signals

[b_LF, a_LF] = fir1(100, 2*[0.5 30]/fs3, "bandpass");

%% Load the envelope reconstruction model weights

useAudChanForDecoding = 7:8;
load(".\models\XX-000\XX000_HighGamma_LowFreq_Model.mat");
gA =  results.NTN_1.gA;

%% Stream and decode!

idx_accu = [];      % This can be commented out later
idx_accu_med = [];  % This can be commented out later
dec_ts_accu = [];   % This can be commented out later

firstBlobObtained = 0;
decodingStarted = 0;

windowSize = 4; %s
hopSize = 0.5; % s
overlap = windowSize - hopSize; % s

noGoodNeuralChannels = 48;
neuralChannels = 9:1:(9 + noGoodNeuralChannels - 1);
neuralChannels_LF = 9:1:(9 + 2 * noGoodNeuralChannels - 1);

%for ctr = 1:1:length(data_accu)
while 1

    [data, ts] = s.get_data();
    %data = data_accu{ctr};
    %ts = ts_accu{ctr};
    
    if firstBlobObtained == 0
        if not(all(ts == 0))
            buffer = data;
            buffer_ts = ts;
            firstBlobObtained = 1;
        end
    else
        
        % Obtain the new samples only since there may be overlap with
        % previous ones

        % Note: the new ts may also contain data in buffer but not with
        % exact time stamps! So don't use intersect function!

        newStartIDX = find(ts > buffer_ts(end), 1);

        if numel(newStartIDX) > 0
            % There is new data
            newData = data(:, newStartIDX:end);
            new_ts  = ts(newStartIDX:end);
            newSamples = numel(new_ts);
        else
            % No new data to add!
            continue;
        end

        if buffer_ts(end) - buffer_ts(1) < 20

            buffer    = [buffer, newData];
            buffer_ts = [buffer_ts, new_ts];

        else
            
            % Get rid of old data and make space for new
            buffer = buffer(:, newSamples+1:end);
            buffer_ts = buffer_ts(newSamples+1:end);

            % Add new data
            buffer    = [buffer, newData];
            buffer_ts = [buffer_ts, new_ts];

        end

        % Perform decoding if you have enough data
        if decodingStarted == 0
            if buffer_ts(end) - buffer_ts(1) >  10 % (windowSize + 0.1) % Changed to a longer duration to keep up with a sharper notch filter
                
                %%%%%%%%%%%%%
                % Do decoding
                %%%%%%%%%%%%%

                % CMR (only on the neural channels)
                buffer_cmr = buffer; % performCMR(buffer, amplifiers);

                % Resample to 400 Hz
                [buffer_cmr_400, buffer_ts_400] = resample(buffer_cmr', buffer_ts, 400);
                buffer_cmr_400 = buffer_cmr_400';

                % HG Env Extraction (only on the neural channels)
                buffer_cmr_400_hg = extractHighGamma(buffer_cmr_400, neuralChannels, b, a, bl, al, bN, aN);

                % LF Signals Extraction (*** NEW ***)
                % buffer_cmr_400_lf = extractLowFrequency(buffer_cmr_400, neuralChannels, b_LF, a_LF);

                % Concat things! (*** NEW ***)
                % buffer_cmr_400_hg_lf = [buffer_cmr_400_hg; buffer_cmr_400_lf(neuralChannels, :)];
                buffer_cmr_400_hg_lf = [buffer_cmr_400_hg; filtfilt(b_LF, a_LF, buffer_cmr_400(neuralChannels, :)')'];
                % buffer_cmr_400_hg_lf = [buffer_cmr_400_hg; buffer_cmr_400(neuralChannels, :)];

                % Resample to 100 Hz
                % [buffer_re, buffer_ts_re] = resample(buffer_cmr_400_hg', buffer_ts_400, 100);
                [buffer_re, buffer_ts_re] = resample(buffer_cmr_400_hg_lf', buffer_ts_400, 100);
                buffer_re = buffer_re';

                % Change neural channels to incorporate LF (*** NEW ***)
                % neuralChannels = 9:1:(9 + 2 * noGoodNeuralChannels - 1);
                
                buffer_re_std = mapstd(buffer_re);  % Normalize each channel
                idx = decodeBuffer(buffer_re_std(:, end - windowSize * 100 + 1:end), neuralChannels_LF, gA, useAudChanForDecoding);

                lastDecodingTimeStep = buffer_ts_re(end) + 0.01;
                idx_accu = [idx_accu, idx];
                idx_accu_med = medfilt1(idx_accu, 3);
                dec_ts_accu = [dec_ts_accu, lastDecodingTimeStep];

                try 
                    curr_idx = idx_accu_med(end-1); % Skip the last one since the one of the elements for computing median is 0.
                catch
                    fprintf("Not enough samples for median filtering! \n");
                    curr_idx = idx_accu_med(end);
                end

                fprintf("\nAAD: %8.2f, IDX: %6.4f, IDX_FILT: %6.4f\n", lastDecodingTimeStep, idx, curr_idx);
                write(u, [lastDecodingTimeStep, curr_idx], "double", "192.168.1.1", 1111);

                % Set decodingStarted flag to 1
                decodingStarted = 1;

            end
        else
            
            if buffer_ts(end) > lastDecodingTimeStep + hopSize

                %%%%%%%%%%%%%
                % Do decoding
                %%%%%%%%%%%%%

                % CMR (only on the neural channels)
                buffer_cmr = buffer; % performCMR(buffer, amplifiers);

                % Resample to 400 Hz
                [buffer_cmr_400, buffer_ts_400] = resample(buffer_cmr', buffer_ts, 400);
                buffer_cmr_400 = buffer_cmr_400';

                % HG Env Extraction (only on the neural channels)
                buffer_cmr_400_hg = extractHighGamma(buffer_cmr_400, neuralChannels, b, a, bl, al, bN, aN);

                % LF Signals Extraction (*** NEW ***)
                % buffer_cmr_400_lf = extractLowFrequency(buffer_cmr_400, neuralChannels, b_LF, a_LF);

                % Concat things! (*** NEW ***)
                % buffer_cmr_400_hg_lf = [buffer_cmr_400_hg; buffer_cmr_400_lf(neuralChannels, :)];
                buffer_cmr_400_hg_lf = [buffer_cmr_400_hg; filtfilt(b_LF, a_LF, buffer_cmr_400(neuralChannels, :)')'];
                % buffer_cmr_400_hg_lf = [buffer_cmr_400_hg; buffer_cmr_400(neuralChannels, :)];

                % Resample to 100 Hz
                % [buffer_re, buffer_ts_re] = resample(buffer_cmr_400_hg', buffer_ts_400, 100);
                [buffer_re, buffer_ts_re] = resample(buffer_cmr_400_hg_lf', buffer_ts_400, 100);
                buffer_re = buffer_re';

                % Change neural channels to incorporate LF (*** NEW ***)
                % neuralChannels = 9:1:(9 + 2 * noGoodNeuralChannels - 1);
                
                buffer_re_std = mapstd(buffer_re);  % Normalize each channel

                %%%%%
                decodeWindowStartIDX = find(buffer_ts_re>=lastDecodingTimeStep-overlap, 1);
                decodeWindowEndIDX = decodeWindowStartIDX + windowSize * 100 - 1;

                idx = decodeBuffer(buffer_re_std(:, decodeWindowStartIDX:decodeWindowEndIDX), neuralChannels_LF, gA, useAudChanForDecoding);

                lastDecodingTimeStep = buffer_ts_re(decodeWindowEndIDX) + 0.01;
                idx_accu = [idx_accu, idx];
                idx_accu_med = medfilt1(idx_accu, 3);
                dec_ts_accu = [dec_ts_accu, lastDecodingTimeStep];

                try 
                    curr_idx = idx_accu_med(end-1); % Skip the last one since the one of the elements for computing median is 0.
                catch
                    fprintf("Not enough samples for median filtering! \n");
                    curr_idx = idx_accu_med(end);
                end

                fprintf("\nAAD: %8.2f, IDX: %6.4f, IDX_FILT: %6.4f\n", lastDecodingTimeStep, idx, curr_idx);
                write(u, [lastDecodingTimeStep, curr_idx], "double", "192.168.1.1", 1111);

            end
            
        end

    end

end

%% Function to save things after running!

save("data_XX_000_" + datestr(datetime, 'mm-dd-yyyy-HH-MM-SS-FFF') + ".mat", 'dec_ts_accu', 'idx_accu', 'idx_accu_med', 'buffer', 'buffer_ts', 'gA');

%% Function to decode

function idx = decodeBuffer(buff, neuralChannels, gA, useAudChanForDecoding)

    lags = -40:0;
    resp = buff(neuralChannels, :);
    [~, rstimA] = StimuliReconstruction([],[],resp',gA,lags);

    stim1 = buff(useAudChanForDecoding(1), :);
    stim2 = buff(useAudChanForDecoding(2), :);

    idx = (corr2(rstimA, stim1) - corr2(rstimA, stim2))/2;

end

%% Function to perform CMR

function newBuffer = performCMR(buffer_in, amplifiers)
    
    newBuffer = buffer_in;

    for ctr = 1:1:length(amplifiers)

        newBuffer(amplifiers{ctr}, :) = newBuffer(amplifiers{ctr}, :) - mean(newBuffer(amplifiers{ctr}, :));

    end

end

%% Function to extract high gamma envelope

function buffer_cmr_400_hg = extractHighGamma(buffer_cmr_400, neuralChannels, b, a, bl, al, bN, aN)

    buffer_cmr_400_hg = buffer_cmr_400;
    
    buffer_cmr_400_hg(neuralChannels, :) = filtfilt(bN, aN, buffer_cmr_400_hg(neuralChannels, :)')';
    buffer_cmr_400_hg(neuralChannels, :) = filtfilt(b, a, buffer_cmr_400_hg(neuralChannels, :)')';
    buffer_cmr_400_hg(neuralChannels, :) = abs( buffer_cmr_400_hg(neuralChannels, :));
    buffer_cmr_400_hg(neuralChannels, :) = filtfilt(bl, al, buffer_cmr_400_hg(neuralChannels, :)')';

end

%% Function to extract low frequency signals

function buffer_cmr_400_lf = extractLowFrequency(buffer_cmr_400, neuralChannels, b_LF, a_LF)

    buffer_cmr_400_lf = buffer_cmr_400;

    buffer_cmr_400_lf(neuralChannels, :) = filtfilt(b_LF, a_LF, buffer_cmr_400_lf(neuralChannels, :)')';

end