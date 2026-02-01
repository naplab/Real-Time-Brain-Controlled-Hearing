%% About

% This scipt runs the OFFLINE phase of the real-time AAD task.

%% Clear the workspace

clc;
close all;
clear all;

%% Change to current directory

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%% Add path

addpath(genpath("./"));

%% Declare subject code

warning("Confirm subject code!");
subjectCode = "XX_000";

%% Declare parameters

fprintf("Please confirm parameters!\n");
TIMESTAMPS = [10, 24, 38];   % Start of Segment 1, 2 and end
sceneSNR = -6; % dB (AAD OFF) case
BG_Difficulty = 15; % [12, 9]; % BG_Difficulty(1) -> Least Difficult -> Low BG Noise -> High supression. % P_Att - P_BG_Noise
MAX_SUPRESSION = 9; % dB - How much to supress un-attend talker in the online phase?
fDecision = 2; % Hz
tDecision = 1/fDecision;

% comparisonQuestionEnabled = 1;

frameLength = 1024;
sampleRate = 44100;

%% Initialise the Markov chain

% IDX = L - R, so +ve => Retain L, Supress R

N_MC = 5;
MC_Neutral = (N_MC + 1)/2; % Neutral state
MC_CurrentState = MC_Neutral; 

Markov_TMRs_2 = linspace(sceneSNR, MAX_SUPRESSION, (N_MC + 1)/2);
Markov_TMRs_1 = linspace(-MAX_SUPRESSION, sceneSNR, (N_MC + 1)/2);

Markov_TMRs = [Markov_TMRs_1, Markov_TMRs_2(2:end)];

alpha_var = 10^(-MAX_SUPRESSION/20);
k_TMRs = 10.^(Markov_TMRs/20);
g2     = sqrt((1 + alpha_var^2)./(1 + k_TMRs.^2));  % Non-Target Gain
g1     = k_TMRs.*g2;                                % Target Gain

% NOTE for case 0, need to determine the gains, such that
% conservation of power is applied.

% Calculate power for 0 dB
g_equal = sqrt((1 + alpha_var^2)./(1 + 1.^2));

%% Load the trial information file

infoTable = readtable("stimuli_design_2024.csv", 'VariableNamingRule', 'preserve');

%% Randomize the trial order

try

    infoTable = readtable(sprintf("logs-csv/%s_start.csv", subjectCode), 'VariableNamingRule', 'preserve');
    fprintf("Table found and loaded for %s.\n", subjectCode);

catch

    % A PTN order has not been generated for this subject. Generate it now.
    fprintf("Generating PTN order for %s...\n", subjectCode);

    % PTN order generation
    PTNOrder = randperm(height(infoTable))';

    % Add to infoTable
    infoTable.PTN = PTNOrder;

    % Sort by PTN, re-arrange columns so that PTN is 2nd column
    infoTable = sortrows(infoTable, "PTN");
    infoTable = infoTable(:, [1, end, 2:end-1]);

    % Write out
    writetable(infoTable, sprintf("logs-csv/%s_start.csv", subjectCode));
    fprintf("PTN order generated.\n");

end

%% Reset list of audio devices

audiodevreset;
info = audiodevinfo;

%% Generate a uifigure where prompting and MCQ happens.

% Make a figure playground
promptFigure = drawUIFigure;

% Prompt cross hair
promptImage(promptFigure, 0);
drawnow;

%% Play trial-by-trial

startTrial = 1; % PTN 
endTrial = 26;  % PTN 

% Print instructions for the experimenter
fprintf("---------------------------------------------\n\n");
fprintf("---------- REAL TIME AAD - OFFLINE ----------\n\n");
fprintf("Please instruct the participant to attend to \nthe speech stream that starts first and ask \nthem to maintain their attention on that \nparticular stream throughout the duration of \nthe trial.\n\n");
fprintf("Please also instruct the participant to click \nthe button upon hearing repeated words in the \nstream that starts first.\n\n");
fprintf("Version: 5.01\n\n");
fprintf("---------------------------------------------\n\n");

x = input("Start the experiment? Input 1 or 0: ");

if x == 1

    % Make a table to store play information
    playInfoTable = cell2table(cell(0, numel(infoTable.Properties.VariableNames) + 4), 'VariableNames', [string(infoTable.Properties.VariableNames), "Start Time", "End Time", "P_BG", "Rating"]);

    for ctr = startTrial:1:endTrial

        % Load trial information
        thisPTN = infoTable{ctr, "PTN"}; assert(thisPTN == ctr);
        thisNTN = infoTable{ctr, "NTN"};
        thisBG = string(infoTable{ctr, "Background"});
        % this_P_BG = input("Enter P_BG: "); %  % Start with -9 dB
        thisBGDifficulty = double(infoTable{ctr, "BG Difficulty"});
        this_P_BG = -BG_Difficulty(thisBGDifficulty);
        thisAttDirn = string(infoTable{ctr, "Attend Direction"});
        thisCondition = infoTable{ctr, "Condition"};

        % Trial conditions
        % -1:  ON -> OFF, System ON, Scene SNR = sceneSNR = -6
        %  0: OFF -> OFF, System OFF, Scene SNR = 0 dB
        % +1: OFF ->  ON, System ON, Scene SNR = sceneSNR = -6

        % rwOnsetTimesAtt = getRWOnsetTimes(thisNTN, 1, rwInfo); % ConvNo = 1 -> Attended stream
        % noRWAttConv = sum(rwOnsetTimesAtt < TIMESTAMPS(end)); % length(fieldnames(rwInfo.("Trial" + string(thisNTN)).Conv1));
        rwOnsetTimesAtt = NaN;
        noRWAttConv = NaN;

        fprintf("PTN %d of %d ...\n", thisPTN, endTrial);
        fprintf("NTN: %d, AttDirn: %s, Condition: %d, P_BG: %d, BG: %s, NoRW: %d\n", thisNTN, thisAttDirn, thisCondition, this_P_BG, thisBG, noRWAttConv);

        pause(2);

        power1 = 0; % dB
        power2 = 0; % dB
        power3 = this_P_BG; % -BG_Difficulty(thisBGDifficulty); %dB;

        scaleFactor1 = 10^(power1/20);
        scaleFactor2 = 10^(power2/20);
        scaleFactor3 = 10^(power3/20) * 1/sqrt(2); % Additional 1/root2 since it is going to be duplicated on both left and right. Want net power = power3.

        % Prompt the cross-hair
        promptFigure.Visible = "on";
        pause(3); 

        % Prompt the direction to attend to with an arrow mark
        if thisAttDirn == "L" % Left

            CURRTIME = 0;
            DIRN = "Left";
            [recHandle2, circHandle, circRadius, circ1Handle, circ2Handle, recBlobHandle1, recBlobHandle2, arrowImgHandle] = drawProgressBar(promptFigure, [0, TIMESTAMPS], CURRTIME, DIRN);
            drawnow;

        else

            CURRTIME = 0;
            DIRN = "Right";
            [recHandle2, circHandle, circRadius, circ1Handle, circ2Handle, recBlobHandle1, recBlobHandle2, arrowImgHandle] = drawProgressBar(promptFigure, [0, TIMESTAMPS], CURRTIME, DIRN);
            drawnow;

        end

        pause(3);

        % Declare the filereaders
        fileReader1 = dsp.AudioFileReader(sprintf("Voice_Mono_UnityRMS_-25dB_2024/Trial_%d_Conv_1.wav", thisNTN), 'SamplesPerFrame', frameLength);
        fileReader2 = dsp.AudioFileReader(sprintf("Voice_Mono_UnityRMS_-25dB_2024/Trial_%d_Conv_2.wav", thisNTN), 'SamplesPerFrame', frameLength);
        fileReader3 = dsp.AudioFileReader(sprintf("BG_Mono_UnityRMS_-25dB/%s.wav", thisBG), 'SamplesPerFrame', frameLength);
        
        % Declare the filewriter (to save what was played)
        fileWriter = dsp.AudioFileWriter(char(sprintf("audio-out/NTN_%d_%s.wav", thisNTN, string(datestr(now,'mm-dd-yyyy_HH-MM-SS')))));

        % Declare the deviceWriter
        deviceWriter = audioDeviceWriter('SampleRate', sampleRate);
        deviceWriter.Device = "MacBook Pro Speakers";
        % deviceWriter.Device = 'External Headphones';
        % deviceWriter.Device = 'Logi Z407';
        % deviceWriter.Device = "Aggregate Device";

        % Reset the Markov state
        if abs(thisCondition) == 1
            MC_CurrentState = MC_Neutral; 
        end

        % Start playing
        
        playStartTime = datetime;
        blocksPlayed = 0;
        lastDecisionTime = 0;
        lastProgBarUpdateTime = 0;
        totalUnderrun = 0;

        while ~(isDone(fileReader1) && isDone(fileReader2))

            trialTime = (blocksPlayed * frameLength)/sampleRate;
            
            % Load segments
            y1 = fileReader1(); % Mono att
            y2 = fileReader2(); % Mono un-att
            y3 = fileReader3(); % Mono BG

            % Apply inital scaling
            y1 = y1*scaleFactor1; % scaleFactor1 = 1 
            y2 = y2*scaleFactor2; % scaleFactor1 = 1 
            y3 = y3*scaleFactor3;


            % Do AAD
            if trialTime - lastDecisionTime > tDecision

                % Make a decision
                if thisCondition == 1 % Off to On

                    if trialTime > TIMESTAMPS(2) % In Segment 2

                        MC_CurrentState = updateMCState(MC_CurrentState, +0.1, N_MC);

                    end

                elseif thisCondition == -1 % On to Off

                    if trialTime > TIMESTAMPS(1) && trialTime < TIMESTAMPS(2) % In Segment 1

                        MC_CurrentState = updateMCState(MC_CurrentState, +0.1, N_MC);

                    elseif trialTime > TIMESTAMPS(2) && MC_CurrentState > MC_Neutral % In Segment 2

                        MC_CurrentState = updateMCState(MC_CurrentState, -0.1, N_MC);

                    end

                end

                % Update the variable
                lastDecisionTime = trialTime;

                % Print the decision
                fprintf("Trial Time: %5.2f, Markov State: %d\n", trialTime, MC_CurrentState);

            end

            % Update the progression bar!
            if trialTime - lastProgBarUpdateTime > 0.2 
                
                updateFrame(recHandle2, circHandle, circRadius, circ1Handle, circ2Handle, recBlobHandle1, recBlobHandle2, trialTime, [0, TIMESTAMPS]);
                drawnow;

                lastProgBarUpdateTime = trialTime;

            end

            % Get the gains
            if abs(thisCondition) == 1

                if strcmp(thisAttDirn, "L")
                    g_L = g1(MC_CurrentState);
                    g_R = g2(MC_CurrentState);
                else
                    g_L = g2(MC_CurrentState);
                    g_R = g1(MC_CurrentState);
                end

            else

                % 0 dB case, with conservation of enery of conversations

                g_L = g_equal; 
                g_R = g_equal; 

            end



            % Play the sound + save the audio out locally also
            switch thisAttDirn

                case "L" % Place the att talker on left

                    numUnderrun = deviceWriter([y1 * g_L + y3, y2 * g_R + y3]); %, y1 * g_L, y2 * g_R]); % [Mix_L, Mix_R, Clean_L, Clean_R]
                                    fileWriter([y1 * g_L + y3, y2 * g_R + y3, y1 * g_L, y2 * g_R]);

                case "R" % Place the att talker on right

                    numUnderrun = deviceWriter([y2 * g_L + y3, y1 * g_R + y3]); %, y2 * g_L, y1 * g_R]); % [Mix_L, Mix_R, Clean_L, Clean_R]
                                    fileWriter([y2 * g_L + y3, y1 * g_R + y3, y2 * g_L, y1 * g_R]);

            end

            blocksPlayed = blocksPlayed + 1;
            totalUnderrun = totalUnderrun + numUnderrun;

        end

        playEndTime = datetime;

        % Update the progression bar to finish
        updateFrame(recHandle2, circHandle, circRadius, circ1Handle, circ2Handle, recBlobHandle1, recBlobHandle2, TIMESTAMPS(end), [0, TIMESTAMPS]);
        drawnow;

        % Print samples underrun
        fprintf("Samples underrun: %d\n", totalUnderrun);
        
        % Release the fileReaders, fileWriters and deviceWriters
        release(fileReader1);
        release(fileReader2);
        release(fileReader3);
        release(fileWriter);
        release(deviceWriter);
        
        % Close the arrow figure
        % close(arrowFigure);

        % Cross hair
        % changeImage(promptFigure, 0); % 0 is crosshair, 1 is left, 2 is rights
        

        % Ask for the rating at the end of the trial
        pause(0);
        if abs(thisCondition) == 1
            rating = promptRating(promptFigure, arrowImgHandle);
        else
            rating = NaN;
        end
        % rating = NaN;

        % Ask for comparison at the end, if PTN is not 1
        % if thisPTN ~= 1 && comparisonQuestionEnabled == 1
        %     compRating = promptComparison(promptFigure);
        % else
        %     compRating = "NaN";
        % end
        % compRating = "NaN";

        % Hide the promptFigure and prepare for the next trial, load cross
        % hair
        promptFigure.Visible = 'off';
        promptImage(promptFigure, 0);
        drawnow;

        % Add to CSV and save the CSV
        toAddRecord = cell2table([table2cell(infoTable(ctr, :)), {datestr(playStartTime, 'mm-dd-yyyy-HH-MM-SS-FFF')}, {datestr(playEndTime, 'mm-dd-yyyy-HH-MM-SS-FFF')}, {this_P_BG}, {rating}], 'VariableNames', [string(infoTable.Properties.VariableNames), "Start Time", "End Time", "P_BG", "Rating"]);
        playInfoTable = [playInfoTable; toAddRecord];
        writetable(playInfoTable, "logs-csv/" + datestr(datetime, 'mm-dd-yyyy-HH-MM-SS-FFF') + ".csv");

        % Pause before the next trial
        if ctr ~= endTrial
                shouldContinue = input("Continue? ");
                if shouldContinue ~= 1
                    break;
                end
        end

    end

end


%% Function to generate a UI figure and return a handle

function promptFigure = drawUIFigure
    
    % Get screen parameters
    screenSize = get(groot, 'ScreenSize');
    screenW = screenSize(3);
    screenH = screenSize(4);
    
    useScreenW = 0.90*screenW;
    useScreenH = 0.90*screenH;
    
    useLeft = 0.05/2*screenW;
    useBottom = 0.125*screenH;
    
    % Generate the figure
    % promptFigure = uifigure('Color', [146, 146, 146]/255, 'Name', 'Prompt', 'Position', [useLeft, useBottom, useScreenW, useScreenH], 'Units', 'pixels', 'Visible', 'off');
    
    % Make a figure
    promptFigure = uifigure('Color', 'White', 'Name', 'Experiment', 'Position', [useLeft, useBottom, useScreenW, useScreenH], 'Units', 'pixels', 'Visible', 'off');
    
end


%% Function to prompt arrows or crosshairs

function promptImage(promptFigure, direction)

    % Clear the ui figure
    clf(promptFigure);

    % Get screen parameters
    screenSize = get(groot, 'ScreenSize');
    screenW = screenSize(3);
    screenH = screenSize(4);
    
    useScreenW = 0.90*screenW;
    useScreenH = 0.90*screenH;
    
    % useLeft = 0.05/2*screenW;
    % useBottom = 0.125*screenH;

    leftArrowURL = "arrow-design/leftArrowWhite.jpeg";
    rightArrowURL = "arrow-design/rightArrowWhite.jpeg";
    crossHairURL = "arrow-design/crossHairWhite.jpeg";

    % Generate the title
    promptTitleLabel = uilabel(promptFigure, "Text", "", "FontSize", 100, "FontWeight", 'bold', 'HorizontalAlignment', 'center', 'Position', [0.15 * useScreenW, 0.85 * useScreenH, 0.7 * useScreenW, 0.2 * useScreenH]);

    % Generate the image
    imgHandle = uiimage(promptFigure, "ScaleMethod", 'fit', 'Position', [0.1 * useScreenW, 0.15 * useScreenH, 0.8 * useScreenW, 0.7 * useScreenH]);

    switch direction

        case 0
            imgHandle.ImageSource = crossHairURL;
        
        case 1 % "L"
            imgHandle.ImageSource = leftArrowURL;

        case 2 % "R"
            imgHandle.ImageSource = rightArrowURL;

    end

end

%% Function to clear elements and draw progress bar

function [recHandle2, circHandle, circRadius, circ1Handle, circ2Handle, recBlobHandle1, recBlobHandle2, arrowImgHandle] = drawProgressBar(figHandle, TIMESTAMPS, CURRTIME, DIRN)

    % Clear the uifigure;
    clf(figHandle);

    % Get the resolution of the figure made
    figPosition = get(figHandle, "Position");
    figWidth = figPosition(3);
    figHeight = figPosition(4);

    % Axes for progression bar
    axWidth = 0.8 * figWidth;
    axHeight = 0.4 * figHeight;
    
    axLeft = 0.1 * figWidth;
    axBottom = 0.5 * figHeight;
    
    % Progression bar properties
    rectBarWidth = TIMESTAMPS(end);
    rectBarHeight = 2/15 * 2.5;
    
    axXLim = [-0.05 * rectBarWidth, 1.05 * rectBarWidth];
    axYLim = [-20, 20];
    
    axHandle = uiaxes(figHandle, 'Position', [axLeft, axBottom, axWidth, axHeight], 'XLim', axXLim, 'DataAspectRatio', [1, 1, 1], 'DataAspectRatioMode', 'manual');
    axHandle.XAxis.Visible = 'off';
    axHandle.YAxis.Visible = 'off';
    axHandle.YLim = [-1, 6];
    
    % Draw full progression bar
    rectHandle = rectangle(axHandle, 'Position', [0, -rectBarHeight/2, rectBarWidth, rectBarHeight]);
    
    % Draw the completion bar
    recHandle2 = rectangle(axHandle, 'Position', [0, -rectBarHeight/2, CURRTIME, rectBarHeight], 'FaceColor', "#77AC30", 'EdgeColor', 'none');

    % Draw lines indicating the major timestamps
    for ctr = 1:1:length(TIMESTAMPS)
    
        tempX = [TIMESTAMPS(ctr), TIMESTAMPS(ctr)];
        tempY = [-rectBarHeight/2*4, +rectBarHeight/2*4];
        line(axHandle, tempX, tempY, 'LineWidth', 4, 'Color', 'k');
    
    end
    
    % Circle for current time
    circRadius = rectBarHeight*1.2;
    circHandle = rectangle(axHandle, 'Position', [CURRTIME - circRadius, -circRadius, 2 * circRadius, 2 * circRadius], 'Curvature', [1, 1], 'FaceColor', '#A2142F');
    
    % Add text to highlight sections 1 and 2
    textHandle1 = text(axHandle, (TIMESTAMPS(2) + TIMESTAMPS(3))/2, rectBarHeight * 15, '1', 'FontSize', 55, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    textHandle2 = text(axHandle, (TIMESTAMPS(3) + TIMESTAMPS(4))/2, rectBarHeight * 15, '2', 'FontSize', 55, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    
    % Circle the text for the current section
    circ1Radius = 1;
    circ1Handle = rectangle(axHandle, 'Position', [(TIMESTAMPS(2) + TIMESTAMPS(3))/2 - circ1Radius, rectBarHeight * 11.7, 2 * circ1Radius, 2 * circ1Radius], 'Curvature', [1, 1], 'FaceColor', 'none', 'LineWidth', 4);
    
    circ2Radius = 1;
    circ2Handle = rectangle(axHandle, 'Position', [(TIMESTAMPS(3) + TIMESTAMPS(4))/2 - circ2Radius, rectBarHeight * 11.7, 2 * circ2Radius, 2 * circ2Radius], 'Curvature', [1, 1], 'FaceColor', 'none', 'LineWidth', 4);

    % Added on Dec. 10 - Blobs to also highlight section number
    recBlobHandle1 = rectangle(axHandle, 'Position', [TIMESTAMPS(2), 3.8, TIMESTAMPS(3) - TIMESTAMPS(2), 2.2], 'Curvature', [0.0, 0.0], 'FaceColor', [1, 0, 0, 0.2], 'EdgeColor', [1, 0, 0, 0]);
    recBlobHandle2 = rectangle(axHandle, 'Position', [TIMESTAMPS(3), 3.8, TIMESTAMPS(4) - TIMESTAMPS(3), 2.2], 'Curvature', [0.0, 0.0], 'FaceColor', [1, 0, 0, 0.2], 'EdgeColor', [1, 0, 0, 0]);
     
    % Control visibility
    if CURRTIME < TIMESTAMPS(2)
        circ1Handle.Visible = 'off';
        circ2Handle.Visible = 'off';
        recBlobHandle1.Visible = 'off';
        recBlobHandle2.Visible = 'off';
    elseif CURRTIME > TIMESTAMPS(2) && CURRTIME < TIMESTAMPS(3)
        circ1Handle.Visible = 'on';
        circ2Handle.Visible = 'off';
        recBlobHandle1.Visible = 'on';
        recBlobHandle2.Visible = 'off';
    elseif CURRTIME < TIMESTAMPS(4)
        circ1Handle.Visible = 'off';
        circ2Handle.Visible = 'on';
        recBlobHandle1.Visible = 'off';
        recBlobHandle2.Visible = 'on';
    end
    
    % Arrow image
    imgWidth = 0.6 * figWidth;
    imgHeight = 0.4 * figHeight;
    
    imgLeft = 0.2 * figWidth;
    imgBottom = 0.1 * figHeight;
    
    arrowImgHandle = uiimage(figHandle, 'Position', [imgLeft, imgBottom, imgWidth, imgHeight], 'ImageSource', "video-frame-design/arrow" + DIRN + ".jpeg");    

end

%% Function to update progression bar

function updateFrame(recHandle2, circHandle, circRadius, circ1Handle, circ2Handle, recBlobHandle1, recBlobHandle2, CURRTIME, TIMESTAMPS)

    % OLD_TIME = CURRTIME - DELTIME;

    if CURRTIME <= TIMESTAMPS(4)
        recHandle2.Position(3) = CURRTIME;
        circHandle.Position(1) = CURRTIME - circRadius;
    end

    % Control visibility
    if CURRTIME < TIMESTAMPS(2)
        circ1Handle.Visible = 'off';
        circ2Handle.Visible = 'off';
        recBlobHandle1.Visible = 'off';
        recBlobHandle2.Visible = 'off';
    elseif CURRTIME > TIMESTAMPS(2) && CURRTIME < TIMESTAMPS(3)
        circ1Handle.Visible = 'on';
        circ2Handle.Visible = 'off';
        recBlobHandle1.Visible = 'on';
        recBlobHandle2.Visible = 'off';
    elseif CURRTIME > TIMESTAMPS(3) && CURRTIME < TIMESTAMPS(4)
        circ1Handle.Visible = 'off';
        circ2Handle.Visible = 'on';
        recBlobHandle1.Visible = 'off';
        recBlobHandle2.Visible = 'on';
    else
        circ1Handle.Visible = 'off';
        circ2Handle.Visible = 'off';
        recBlobHandle1.Visible = 'off';
        recBlobHandle2.Visible = 'off';
    end

end

%% Function to change an image

function changeImage(promptFigure, direction)

    leftArrowURL = "arrow-design/leftArrowWhite.jpeg";
    rightArrowURL = "arrow-design/rightArrowWhite.jpeg";
    crossHairURL = "arrow-design/crossHairWhite.jpeg";

    switch direction

        case 0
            targetURL = crossHairURL;
        
        case 1 % "L"
            targetURL = leftArrowURL;

        case 2 % "R"
            targetURL = rightArrowURL;

    end

    for ctr = 1:1:length(promptFigure.Children)

        if strcmp(promptFigure.Children(ctr).Type, "uiimage")
            promptFigure.Children(ctr).ImageSource = targetURL;
        end

    end

end

%% Function to prompt rating question

function rating = promptRating(promptFigure, arrowImgHandle)

    % Get the position of arrowImgHandle
    thisPosition = arrowImgHandle.Position;

    % Delete the arrowImgHandle
    delete(arrowImgHandle);

    % Generate question label
    ratingQuestionLabel = uilabel(promptFigure, "Text", "Which part was more DIFFICULT to follow?", "FontSize", 38, "FontWeight", 'bold', 'HorizontalAlignment', 'center', 'Position', [thisPosition(1), thisPosition(2) + thisPosition(4)/2, thisPosition(3), thisPosition(4)/2]);

    % Generate button group
    % bgWidth = 0.95 * useScreenW;
    % bgHeight = 0.2 * useScreenH;

    bg = uibuttongroup(promptFigure, 'Position', [thisPosition(1), thisPosition(2)/3, thisPosition(3), thisPosition(4) * 0.75], 'SelectionChangedFcn', @selectionMade, 'BackgroundColor', [1, 1, 1], 'BorderType', 'none');

    % Generate buttons

    bFontSize = 64;

    b1 = uiradiobutton(bg, "Text", " 1", "Position", [1 * thisPosition(3)/3 - thisPosition(3)/20, thisPosition(4)/2, thisPosition(3)/4, thisPosition(4)/4], 'FontSize', bFontSize, 'FontWeight', 'bold');
    b2 = uiradiobutton(bg, "Text", " 2", "Position", [2 * thisPosition(3)/3 - thisPosition(3)/20, thisPosition(4)/2, thisPosition(3)/4, thisPosition(4)/4], 'FontSize', bFontSize, 'FontWeight', 'bold');
    b3 = uiradiobutton(bg, "Text", " 3", "Position", [3 * thisPosition(3)/3 - thisPosition(3)/20, thisPosition(4)/2, thisPosition(3)/4, thisPosition(4)/4], 'FontSize', bFontSize, 'FontWeight', 'bold', 'Value', 1, 'Visible', 'off');

    % leftShift = 0.075/2 * bgWidth;
    % bottomShift = 0.125 * bgHeight;
    % buttonWidth = 0.18 * bgWidth;
    % buttonHeigth = 0.75 * bgHeight;
    % 
    % shiftBetweenButtons = 0.2 * bgWidth;

    % fprintf("I am here!\n");

    % Wait for response
    uiwait(promptFigure);
    
    % Capture the rating
    ratingText = bg.SelectedObject.Text;
    rating = double(string(ratingText(2)));
    
    % Close the figure
    % close(promptFigure);

end

%% Function to prompt comparison rating

function compRating = promptComparison(promptFigure)

    % Clear the ui figure
    clf(promptFigure);

    % Get screen parameters
    screenSize = get(groot, 'ScreenSize');
    screenW = screenSize(3);
    screenH = screenSize(4);
    
    useScreenW = 0.95*screenW;
    useScreenH = 0.75*screenH;
    
    useLeft = 0.05/2*screenW;
    useBottom = 0.125*screenH;

    % Generate figure
    % ratingFigure = uifigure('Color', 'white', 'Name', 'Rating Prompt', 'Position', [useLeft, useBottom, useScreenW, useScreenH], 'Units', 'pixels');

    % Generate label
    ratingQuestionLabel = uilabel(promptFigure, "Text", "Compared to the previous trial, how difficult was the current trial?", "FontSize", 38, "FontWeight", 'bold', 'HorizontalAlignment', 'center', 'Position', [0.1 * useScreenW, 0.7 * useScreenH, 0.8 * useScreenW, 0.2 * useScreenH]);

    % Generate button group
    bgWidth = 0.8 * useScreenW;
    bgHeight = 0.2 * useScreenH;

    bg = uibuttongroup(promptFigure, 'Position', [0.1 * useScreenW, 0.5 * useScreenH, bgWidth, bgHeight], 'SelectionChangedFcn', @selectionMade, 'BackgroundColor', [146, 146, 146]/255);

    % Generate buttons

    leftShift = 0.19 * bgWidth;
    bottomShift = 0.125 * bgHeight;
    buttonWidth = 0.13 * bgWidth;
    buttonHeigth = 0.75 * bgHeight;

    shiftBetweenButtons = 0.25 * bgWidth;

    bFontSize = 32;

    b1 = uiradiobutton(bg, "Text", " Harder", "Position", [leftShift, bottomShift, buttonWidth, buttonHeigth], 'FontSize', bFontSize, 'FontWeight', 'bold');
    b2 = uiradiobutton(bg, "Text", " Same", "Position", [shiftBetweenButtons + leftShift, bottomShift, buttonWidth, buttonHeigth], 'FontSize', bFontSize, 'FontWeight', 'bold');
    b3 = uiradiobutton(bg, "Text", " Easier", "Position", [2*shiftBetweenButtons+leftShift, bottomShift, buttonWidth, buttonHeigth], 'FontSize', bFontSize, 'FontWeight', 'bold');
    % b4 = uiradiobutton(bg, "Text", " 4. Easy", "Position", [3*shiftBetweenButtons+leftShift, bottomShift, buttonWidth, buttonHeigth], 'FontSize', bFontSize, 'FontWeight', 'bold');
    % b5 = uiradiobutton(bg, "Text", " 5. Very Easy", "Position", [4*shiftBetweenButtons+leftShift, bottomShift, buttonWidth, buttonHeigth], 'FontSize', bFontSize, 'FontWeight', 'bold');
    b6 = uiradiobutton(bg, "Text", " Skip", "Position", [3.75*shiftBetweenButtons, bottomShift/10, buttonWidth/4, buttonHeigth/4], 'FontSize', 10, 'FontWeight', 'bold');
    b6 = uiradiobutton(bg, "Text", " Initial", "Position", [3*shiftBetweenButtons, 0.775 * bgHeight, buttonWidth/4, buttonHeigth/4], 'FontSize', 10, 'FontWeight', 'bold', 'Value', 1, 'Visible', 'off');

    % fprintf("I am here!\n");

    % Wait for response
    uiwait(promptFigure);
    
    % Capture the rating
    ratingText = bg.SelectedObject.Text;
    compRating = string(ratingText(2));

    if ~ismember(compRating, ["H", "S", "E"])
        compRating = "NaN";
    end
    
    % Close the figure
    % close(ratingFigure);

end

%% Function to resume execution once a selection has been made

function selectionMade(src, event)

    uiresume(src.Parent);

end

%% Funtion to get RW Onset Times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions to get RW onset times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rwOnsetTimes = getRWOnsetTimes(NTN, convNo, rwInfo)

    noRWs = length(fieldnames(rwInfo.("Trial" + string(NTN)).("Conv" + string(convNo))));
    
    rwOnsetTimes = zeros(1, noRWs);
    
    for ctr = 1:1:noRWs
        
        rwOnsetTimes(ctr) = rwInfo.("Trial" + string(NTN)).("Conv" + string(convNo)).("RW" + string(ctr)).onset;
        
    end
    
end

%% Function to update MC state

function MC_CurrentState = updateMCState(MC_CurrentState, IDX, MC_N_Levels)

    if IDX > 0
        currState = MC_CurrentState + 1;
    elseif IDX < 0
        currState = MC_CurrentState - 1;
    else
        currState = MC_CurrentState;
    end

    if currState > MC_N_Levels
        MC_CurrentState = MC_N_Levels;
    elseif currState < 1
        MC_CurrentState = 1;
    else
        MC_CurrentState = currState;
    end

end