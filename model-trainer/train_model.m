%% About

% Train brain -> attended conversation envelope reconstruction models 
% LF (low-freq)
% HG (high-gamma)
% LF + HG

%% Prepare the workspace

clc;
clear all;

%% Change to current directory and add path to preprocessing functions

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath("ECoGDataPipeline/"));

% General functions
addpath(genpath("./functions/"));

%% Declare subjects

thisSubjectPLs = ["XX"];
thisSubjectIDs = ["000"];

subjectCTR = 1;

thisSubjectPL = thisSubjectPLs(subjectCTR);
thisSubjectID = thisSubjectIDs(subjectCTR);

%% Load the htkraw out struct - out3 (400 Hz)

dpath = ['../neural-data/' char(thisSubjectPL) '-' char(thisSubjectID) filesep];
readFileName = sprintf("out_%s-%s_B1_htkraw.mat", thisSubjectPL, thisSubjectID);
load([dpath filesep 'processed/B1/' char(readFileName)], 'out');

out3 = out;
clear out;

%% Rearrange based on NTN

NTNs_accu = [];

for ctr = 1:1:length(out3)
    
    thisNTN = string(out3(ctr).name);
    thisNTN = split(thisNTN, "_");
    thisNTN = double(thisNTN(2));
    
    NTNs_accu = [NTNs_accu, thisNTN];
    
end

[temp, idx] = sort(NTNs_accu);
out3 = out3(idx);

%% Declare parameters

fNeuralData = 100; % Hz % Target sampling rate

%% Load the speech responsive electrodes

load([dpath filesep 'processed/B1/responsive_electrodes.mat']);

%% Extract high-gamma for out3

fs3 = 400; % Hz
notchFreq = 120; % Hz
[b,a] = fir1(100,2*[70 150]/fs3);
[bl, al] = fir1(100, 2*20/fs3); % 20 Hz LPF
[bN, aN] = butter(5, [notchFreq - 1 notchFreq + 1]/(fs3/2), "stop");


for ctr = 1:1:length(out3)

    temp = out3(ctr).resp;
    temp = filtfilt(bN, aN, temp')';
    temp = filtfilt(b, a, temp')';
    temp = abs(temp);
    temp = filtfilt(bl, al, temp')';
    temp = resample(temp', fNeuralData, fs3)'; % Resample down to 100 Hz
    out3(ctr).resp_hg = temp;    
    out3(ctr).dataf_hg = fNeuralData;

end

%% Extract low-frequency for out3

[b_LF, a_LF] = fir1(100, 2*[0.5 30]/fs3, "bandpass");

for ctr = 1:1:length(out3)

    temp = out3(ctr).resp;
    temp = filtfilt(b_LF, a_LF, temp')';
    % temp = filtfilt(bpFilter, temp')';
    temp = resample(temp', fNeuralData, fs3)'; % Resample down to 100 Hz
    out3(ctr).resp_lf = temp; 
    out3(ctr).dataf_lf = fNeuralData;


end

%% Classifier details preparation

% Load attend directions
stimTable = readtable(sprintf("../stimuli-design/offline/logs-csv/%s_%s_end.csv", thisSubjectPL, thisSubjectID), 'VariableNamingRule', 'preserve');
attDirn = []; % Sorted by PTN, needs to be sorted by NTN
NTNAccu = [];
condnAccu = [];

fprintf("CHECK CHANNEL CONVENTION!\n");

for ctr = 1:1:size(stimTable, 1)
    
    if strcmp(string(stimTable{ctr, "Attend Direction"}), "L")
        attDirn = [attDirn, 1]; % Left channel was connected to channel 1
    else
        attDirn = [attDirn, 2]; % Right channel was connected to channel 2
    end

    NTNAccu = [NTNAccu, stimTable{ctr, "NTN"}];
    condnAccu = [condnAccu, stimTable{ctr, "Condition"}];
    
end

% Sorting by NTN
[NTNAccu, sortIDX] = sort(NTNAccu);
attDirn = attDirn(sortIDX);
condnAccu = condnAccu(sortIDX);

silenceDurnBeforeStimulusOnset = 0.5;               % [in seconds]
skipDurnPostStimuliOnset = 4;                       % [in seconds] NOTE: 3 s is ST, additional 1 s to allow for settling post MT introduction
silenceDurnAfterStimulusOffset = 0.5;               % [in seconds] % If this value = 1.5, then 0.5 s = post stim silence, 1.0 s = etches part of the trial.

fNeuralData = 100;                                  % Hz sampling rate of neural data

skipStartSamples = (silenceDurnBeforeStimulusOnset + skipDurnPostStimuliOnset) * fNeuralData;
skipEndSamples = (silenceDurnAfterStimulusOffset) * fNeuralData;

%% Design filters for waveform envelope extraction - FOR RIPPLE

fsAnalog = 1000;

% Desing a filter
h  = fdesign.lowpass('N,F3dB', 4, 8, fsAnalog); % 4th order low pass filter with a 3 dB cut-off at 8 Hz
hObj = design(h, 'butter');

%% Prepare the data, dividing 1 and 2 segments!

toPassNeuralResponse = {};
toPassAttEnv = {};
toPassUnAttEnv = {};
toPassNTNs = [];

neuralCombo_LF = [1, 0, 1];
neuralCombo_HF = [0, 1, 1];

neuralComboCTR = 3;

USE_LF = neuralCombo_LF(neuralComboCTR);
USE_HG = neuralCombo_HF(neuralComboCTR);

% CHANGE THIS TO 1:1:... in the final version 
fprintf("Verify training trials range! \n");

usefulTrials = 1:1:26;

warning("Removing bad trials for this subject!");
usefulTrials([14, 23]) = []; % Remove bad trials

toPassTrials = [];
orgNTNaccu = [];
passTrialCTR = 0;

for NTN = usefulTrials
    
    thisNTN = string(out3(NTN).name);
    thisNTN = split(thisNTN, "_");
    thisNTN = double(thisNTN(2));
    
    assert(thisNTN == NTN);

    % Process the envelopes 

    % TDT Env

    % tempAttEnv = out3(NTN).audrecord(:, attDirn(NTN) + 14);
    % tempUnAttEnv = out3(NTN).audrecord(:, 3 - attDirn(NTN) + 14);
    % 
    % fsEnv = 50*10^6/2^13;
    % fprintf("\n Check fs of neural data and audio!");
    % [p, q] = rat(100/fsEnv, 1e-8);
    % 
    % tempAttEnv = resample(tempAttEnv, p, q)';
    % tempUnAttEnv = resample(tempUnAttEnv, p, q)';

    % Ripple Env

    warning("Modified for Ripple!");

    tempAttEnv = getAudioEnvRipple(out3(NTN).audrecord(:, attDirn(NTN) + 2), hObj);
    tempUnAttEnv = getAudioEnvRipple(out3(NTN).audrecord(:, 3 - attDirn(NTN) + 2), hObj);

    fsEnv = 1000; % 50*10^6/2^13;
    
    [p, q] = rat(100/fsEnv, 1e-8);

    tempAttEnv = resample(tempAttEnv, p, q)';
    tempUnAttEnv = resample(tempUnAttEnv, p, q)';


    % MATLAB Env

    % tempAttEnv = getAudioEnvelope(out3(NTN).sound(:, attDirn(NTN)), out3(NTN).soundf);
    % tempAttEnv  = resample(tempAttEnv, fNeuralData, out3(NTN).soundf)';
    % 
    % tempUnAttEnv = getAudioEnvelope(out3(NTN).sound(:, 3 - attDirn(NTN)), out3(NTN).soundf);
    % tempUnAttEnv  = resample(tempUnAttEnv, fNeuralData, out3(NTN).soundf)';

    % MATLAB Env - HF, LF

    % tempAttEnv = [getAudioEnvelope_LF(out3(NTN).sound(:, attDirn(NTN)), out3(NTN).soundf), getAudioEnvelope_HF(out3(NTN).sound(:, attDirn(NTN)), out3(NTN).soundf)]'; % 2 x T
    % tempAttEnv = resample(tempAttEnv', fNeuralData, out3(NTN).soundf)';
    % 
    % tempUnAttEnv = [getAudioEnvelope_LF(out3(NTN).sound(:, 3 - attDirn(NTN)), out3(NTN).soundf), getAudioEnvelope_HF(out3(NTN).sound(:, 3 - attDirn(NTN)), out3(NTN).soundf)]'; % 2 x T
    % tempUnAttEnv = resample(tempUnAttEnv', fNeuralData, out3(NTN).soundf)';

    % Spectrogram

    % tempAttEnv = getAuditorySpectrogramMagnitude(out3(NTN).sound(:, attDirn(NTN)), out3(NTN).soundf, 0, 0);
    % tempAttEnv = tempAttEnv.^0.05;
    % 
    % tempUnAttEnv = getAuditorySpectrogramMagnitude(out3(NTN).sound(:, 3 - attDirn(NTN)), out3(NTN).soundf, 0, 0);
    % tempUnAttEnv = tempUnAttEnv.^0.05;

    try 
        tempAttEnv = tempAttEnv(:, 1:size(out3(NTN).resp_hg, 2));
        tempUnAttEnv = tempUnAttEnv(:, 1:size(out3(NTN).resp_hg, 2));
    catch
        warning("Stimulus shorter than neural response!\n");

        toPad = size(out3(NTN).resp_hg, 2) - size(tempAttEnv, 2);
        tempAttEnv = [tempAttEnv, zeros(size(tempAttEnv, 1), toPad)];
        tempUnAttEnv = [tempUnAttEnv, zeros(size(tempUnAttEnv, 1), toPad)];

    end

    % Generate time stamps
    t = generateTimePoints(out3(1).resp_hg, 100, 0);

    s1 = find(t > 10, 1);
    e1 = find(t >= 24, 1);

    s2 = find(t > 24, 1);
    e2 = find(t >= 38, 1);

    % Switch things based on label, FLIP FOR AAD OFF


    switch condnAccu(NTN)

        case 1 % AAD Off to On, FLIP FOR AAD OFF => SEG 1

            if USE_LF == 0 && USE_HG == 1
            
                seg1NeuralResponse = out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s1:e1);
                seg2NeuralResponse = out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s2:e2);
            
            elseif USE_LF == 1 && USE_HG == 0
            
                seg1NeuralResponse = out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s1:e1);
                seg2NeuralResponse = out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s2:e2);
            
            else % USE_LF == 1 && USE_HG == 1
            
                seg1NeuralResponse = [out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s1:e1); out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s1:e1)];
                seg2NeuralResponse = [out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s2:e2); out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s2:e2)];
            
            end

            seg1AttEnv = tempAttEnv(:, s1:e1);
            seg2AttEnv = tempAttEnv(:, s2:e2);

            seg1UnAttEnv = tempUnAttEnv(:, s1:e1);
            seg2UnAttEnv = tempUnAttEnv(:, s2:e2);

            % Pass data - Seg 1
            toPassNeuralResponse{end+1} = seg1NeuralResponse;
            toPassAttEnv{end+1} = seg1AttEnv; % seg1UnAttEnv; % FLIP
            toPassUnAttEnv{end+1} = seg1UnAttEnv; % seg1AttEnv; % FLIP
            passTrialCTR = passTrialCTR + 1;
            toPassTrials = [toPassTrials, passTrialCTR];
            orgNTNaccu = [orgNTNaccu, NTN];

            % Pass data - Seg 2
            toPassNeuralResponse{end+1} = seg2NeuralResponse;
            toPassAttEnv{end+1} = seg2AttEnv;
            toPassUnAttEnv{end+1} = seg2UnAttEnv;
            passTrialCTR = passTrialCTR + 1;
            toPassTrials = [toPassTrials, passTrialCTR];
            orgNTNaccu = [orgNTNaccu, NTN];

        case -1 % AAD On to Off, FLIP FOR AAD OFF => SEG 2

            if USE_LF == 0 && USE_HG == 1
            
                seg1NeuralResponse = out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s1:e1);
                seg2NeuralResponse = out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s2:e2);
            
            elseif USE_LF == 1 && USE_HG == 0
            
                seg1NeuralResponse = out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s1:e1);
                seg2NeuralResponse = out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s2:e2);
            
            else % USE_LF == 1 && USE_HG == 1
            
                seg1NeuralResponse = [out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s1:e1); out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s1:e1)];
                seg2NeuralResponse = [out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, s2:e2); out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, s2:e2)];
            
            end

            seg1AttEnv = tempAttEnv(:, s1:e1);
            seg2AttEnv = tempAttEnv(:, s2:e2);

            seg1UnAttEnv = tempUnAttEnv(:, s1:e1);
            seg2UnAttEnv = tempUnAttEnv(:, s2:e2);

            % Pass data - Seg 1
            toPassNeuralResponse{end+1} = seg1NeuralResponse;
            toPassAttEnv{end+1} = seg1AttEnv; 
            toPassUnAttEnv{end+1} = seg1UnAttEnv; 
            passTrialCTR = passTrialCTR + 1;
            toPassTrials = [toPassTrials, passTrialCTR];
            orgNTNaccu = [orgNTNaccu, NTN];

            % Pass data - Seg 2
            toPassNeuralResponse{end+1} = seg2NeuralResponse;
            toPassAttEnv{end+1} = seg2AttEnv; % seg2UnAttEnv; % FLIP
            toPassUnAttEnv{end+1} = seg2UnAttEnv; % seg2AttEnv; % FLIP
            passTrialCTR = passTrialCTR + 1;
            toPassTrials = [toPassTrials, passTrialCTR];
            orgNTNaccu = [orgNTNaccu, NTN];

        case 0 % 0 dB SNR -> Fully flip to unattended

            % toPassNeuralResponse{end+1} = out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, skipStartSamples+1:end-skipEndSamples);
            if USE_LF == 0 && USE_HG == 1
            
                toPassNeuralResponse{end+1} = out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, skipStartSamples+1:end-skipEndSamples);
            
            elseif USE_LF == 1 && USE_HG == 0

                toPassNeuralResponse{end+1} = out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, skipStartSamples+1:end-skipEndSamples);
            
            else % USE_LF == 1 && USE_HG == 1
            
                toPassNeuralResponse{end+1} = [out3(NTN).resp_hg(thisSubjectResponsiveElectrodes, skipStartSamples+1:end-skipEndSamples); out3(NTN).resp_lf(thisSubjectResponsiveElectrodes, skipStartSamples+1:end-skipEndSamples)];

            end

            toPassAttEnv{end+1} = tempAttEnv(:, skipStartSamples+1:end-skipEndSamples); % FLIP 
            toPassUnAttEnv{end+1} = tempUnAttEnv(:, skipStartSamples+1:end-skipEndSamples); % FLIP 

            % toPassAttEnv{end+1} = tempUnAttEnv(skipStartSamples+1:end-skipEndSamples); % FLIP 
            % toPassUnAttEnv{end+1} = tempAttEnv(skipStartSamples+1:end-skipEndSamples); % FLIP 

            passTrialCTR = passTrialCTR + 1;
            toPassTrials = [toPassTrials, passTrialCTR];
            orgNTNaccu = [orgNTNaccu, NTN];

    end
    
    
end

toPassNeuralResponse = normaliseData(toPassNeuralResponse);
toPassAttEnv = cellfun(@mapstd, toPassAttEnv, 'UniformOutput', false);
toPassUnAttEnv = cellfun(@mapstd, toPassUnAttEnv, 'UniformOutput', false);

% toPassAttEnv = normaliseData(toPassAttEnv);
% toPassUnAttEnv = normaliseData(toPassUnAttEnv);

%% Train decoders, gA gU

if USE_LF == 0 && USE_HG == 1

    resultsFileName = sprintf("%s-%d-HG", thisSubjectPL, thisSubjectID);

elseif USE_LF == 1 && USE_HG == 0

    resultsFileName = sprintf("%s-%d-LF", thisSubjectPL, thisSubjectID);

else

    resultsFileName = sprintf("%s-%d-LF-HG", thisSubjectPL, thisSubjectID);

end

for NTN = toPassTrials
    
    fprintf("LO TO - NTN %d...\n", NTN);
    
    testNTN = NTN;
    testNTN_IDX = find(toPassTrials == testNTN);
    
    trainNTNs = setdiff(toPassTrials, testNTN);
    
    lags = [-40:0];
    resp = []; stimA = []; stimU = [];
    
    for ctr = trainNTNs
        
        trainNTN_IDX = find(toPassTrials == ctr);
        
        tmp = toPassNeuralResponse{trainNTN_IDX};
        resp = [resp tmp zeros(size(tmp,1),.1*fNeuralData)];

        tmp = toPassAttEnv{trainNTN_IDX};
        stimA = [stimA tmp zeros(size(tmp,1),.1*fNeuralData)];

        tmp = toPassUnAttEnv{trainNTN_IDX};
        stimU = [stimU tmp zeros(size(tmp,1),.1*fNeuralData)];
        
    end
    
    [this_gA, ~] = StimuliReconstruction(stimA,resp',[],[],lags);
    [this_gU, ~] = StimuliReconstruction(stimU,resp',[],[],lags);
    results.("NTN_" + string(NTN)).gA = this_gA;
    results.("NTN_" + string(NTN)).gU = this_gU;
    
    % Run on test NTN
    [~, rstimA] = StimuliReconstruction([],[],toPassNeuralResponse{testNTN_IDX}',this_gA,lags);
    [~, rstimU] = StimuliReconstruction([],[],toPassNeuralResponse{testNTN_IDX}',this_gU,lags);
    
    results.("NTN_" + string(NTN)).rstimA = rstimA;
    results.("NTN_" + string(NTN)).rstimU = rstimU;
    results.("NTN_" + string(NTN)).stimA = toPassAttEnv{testNTN_IDX};
    results.("NTN_" + string(NTN)).stimU = toPassUnAttEnv{testNTN_IDX};
    
end

% Save results
cust_mkdir(sprintf("results/%s-%d/", thisSubjectPL, thisSubjectID));
save(sprintf("results/%s-%d/%s.mat", thisSubjectPL, thisSubjectID, resultsFileName), 'results');

load handel
sound(y,Fs);

%% Accumulate results as a function of window size

winSizes = 0.5:0.5:12;
overlap = 50;

metric1_accu = NaN(1, length(winSizes));
metric2_accu = NaN(1, length(winSizes));

for winSizeCTR = 1:1:length(winSizes)
    
    thisWindowSize = winSizes(winSizeCTR);
    
    [corrAttEst_Att, corrAttEst_UnAtt, corrUnAttEst_Att, corrUnAttEst_UnAtt, windowLabels] = getCorrelations_gA_gU(results, thisWindowSize, overlap, fNeuralData);

    % Get Accuracies

    metric1Accu = sum(corrAttEst_Att > corrAttEst_UnAtt)/length(corrAttEst_Att) * 100; % gA

    [correctPred, incorrectPred] = metric2(corrAttEst_Att(:, 1), corrAttEst_UnAtt(:, 1), corrUnAttEst_Att(:, 1), corrUnAttEst_UnAtt(:, 1));
    metric2Accu = sum(correctPred>incorrectPred)/numel(correctPred)*100; % gA, gU

    metric1_accu(winSizeCTR) = metric1Accu;
    metric2_accu(winSizeCTR) = metric2Accu;

    fprintf("Window Size: %.2f s, Tot Windows = %d, gA Accuracy: %.1f, gA gU Accuracy: %.1f \n", thisWindowSize, length(corrAttEst_Att), metric1Accu, metric2Accu);
    
end

%% Plot the stuff

figure; 
plot(winSizes, metric1_accu, '-*', 'LineWidth', 2, 'DisplayName', 'gA'); hold on;
plot(winSizes, metric2_accu, '-*', 'LineWidth', 2, 'DisplayName', 'gA + gU');
xlabel("Window Size (s)");
ylabel("Decoding Accuracy (%)");
legend();
set(gca, 'FontSize', 20);
grid on;


%% Save decoding plot data - High Gamma

decodingPlotData = struct;

decodingPlotData.HG.winSizes = winSizes;
decodingPlotData.HG.metric1_accu = metric1_accu;
decodingPlotData.HG.metric2_accu = metric2_accu;

%% Save decoding plot data - High Gamma + LF

decodingPlotData.HG_LF.winSizes = winSizes;
decodingPlotData.HG_LF.metric1_accu = metric1_accu;
decodingPlotData.HG_LF.metric2_accu = metric2_accu;

%% Function that returns correlations for all windows

function [corrAttEst_Att, corrAttEst_UnAtt, windowLabels] = getCorrelations(results, thisWindowSize, overlap, fNeuralData)

    corrAttEst_Att = [];
    corrAttEst_UnAtt = [];
    windowLabels = [];
    
    % Processing preparation
    hopSize = (1 - overlap/100) * thisWindowSize * fNeuralData;
    
    availNTNNames = string(fieldnames(results));
    
    for ctr = 1:1:length(availNTNNames)
        
        thisNTNName = availNTNNames(ctr);
        
        stimA = results.(thisNTNName).stimA;
        stimU = results.(thisNTNName).stimU;
        
        rstimA = results.(thisNTNName).rstimA;
        
        windowStartIDXs = 1:hopSize:numel(stimA);
        
        winCTR = 1;
        
        for windowStartIDX = windowStartIDXs
            
            windowEndIDX = windowStartIDX + thisWindowSize * fNeuralData - 1;
            
            if windowEndIDX < numel(stimA)
                corrAttEst_Att = [corrAttEst_Att; corr2(stimA(windowStartIDX:windowEndIDX), rstimA(windowStartIDX:windowEndIDX))];
                corrAttEst_UnAtt = [corrAttEst_UnAtt; corr2(stimU(windowStartIDX:windowEndIDX), rstimA(windowStartIDX:windowEndIDX))];
                windowLabels = [windowLabels; thisNTNName + "_" + string(winCTR)];
                winCTR = winCTR + 1;
            end
            
        end
        
    end

end

%% Function that returns correlations for all windows (support for gA, gU)

function [corrAttEst_Att, corrAttEst_UnAtt, corrUnAttEst_Att, corrUnAttEst_UnAtt, windowLabels] = getCorrelations_gA_gU(results, thisWindowSize, overlap, fNeuralData)

    corrAttEst_Att = [];
    corrAttEst_UnAtt = [];

    corrUnAttEst_Att = [];
    corrUnAttEst_UnAtt = [];

    windowLabels = [];
    
    % Processing preparation
    hopSize = (1 - overlap/100) * thisWindowSize * fNeuralData;
    
    availNTNNames = string(fieldnames(results));
    
    for ctr = 1:1:length(availNTNNames)
        
        thisNTNName = availNTNNames(ctr);
        
        stimA = results.(thisNTNName).stimA;
        stimU = results.(thisNTNName).stimU;
        
        rstimA = results.(thisNTNName).rstimA;
        rstimU = results.(thisNTNName).rstimU;
        
        windowStartIDXs = 1:hopSize:size(stimA, 2);
        
        winCTR = 1;
        
        for windowStartIDX = windowStartIDXs
            
            windowEndIDX = windowStartIDX + thisWindowSize * fNeuralData - 1;
            
            if windowEndIDX < size(stimA, 2)
                corrAttEst_Att = [corrAttEst_Att; corr2(stimA(:, windowStartIDX:windowEndIDX), rstimA(:, windowStartIDX:windowEndIDX))];
                corrAttEst_UnAtt = [corrAttEst_UnAtt; corr2(stimU(:, windowStartIDX:windowEndIDX), rstimA(:, windowStartIDX:windowEndIDX))];

                corrUnAttEst_Att = [corrUnAttEst_Att; corr2(stimA(:, windowStartIDX:windowEndIDX), rstimU(:, windowStartIDX:windowEndIDX))];
                corrUnAttEst_UnAtt = [corrUnAttEst_UnAtt; corr2(stimU(:, windowStartIDX:windowEndIDX), rstimU(:, windowStartIDX:windowEndIDX))];


                windowLabels = [windowLabels; thisNTNName + "_" + string(winCTR)];
                winCTR = winCTR + 1;
            end
            
        end
        
    end

end

%% gA gU metric

% Metric 2, Nima's

function [correctPred, incorrectPred] = metric2(corrAttEst_Att, corrAttEst_UnAtt, corrUnAttEst_Att, corrUnAttEst_UnAtt)

    correctPred = corrAttEst_Att - corrUnAttEst_Att;
    incorrectPred = corrAttEst_UnAtt - corrUnAttEst_UnAtt;

end

%% Audio Envelope - LF

function env = getAudioEnvelope_LF(y, fs)

    % BPF
    [b, a] = fir1(100, 4000/(fs/2), "low");

    % LPF
    [bl, al] = fir1(100, 10/(fs/2), "low"); % 10 Hz LPF

    y_bpf = filtfilt(b, a, y);
    y_bpf = abs(y_bpf);
    env = filtfilt(bl, al, y_bpf);

end

%% Audio Envelope - HF

function env = getAudioEnvelope_HF(y, fs)

    % BPF
    [b, a] = fir1(100, [4000, 7900]/(fs/2), "bandpass");

    % LPF
    [bl, al] = fir1(100, 10/(fs/2), "low"); % 10 Hz LPF

    y_bpf = filtfilt(b, a, y);
    y_bpf = abs(y_bpf);
    env = filtfilt(bl, al, y_bpf);

end

%% Function to get audio envelope

% Extract the envelope and apply the filter

function audEnv = getAudioEnvRipple(inData, hObj)

    % inData: (channels, time)

    % Extract the envelope
    y_hilb_env = abs(hilbert(inData')); % (time, channels)

    % Apply the filter
    audEnv = filtfilthd(hObj, y_hilb_env); % (time, channels)

    % Transpose
    audEnv = audEnv'; % (channels, time)


end

