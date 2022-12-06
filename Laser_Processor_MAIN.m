%% LA-ICPMS data processing routine developed by Nathaniel Mollica
% Contact info:
% Woods Hole Oceanographic Institution
% 266 Woods Hole Road MS#23
% Woods Hole, MA 02360
% nmollica@whoi.edu

% Version: 1.0
% Last modified: 12/5/2022

clear; close all force; clc

global UI tabgp tab folder

PrimaryStandardValues.BCa  = 459.6e-6;
PrimaryStandardValues.MgCa = 4.199e-3;
PrimaryStandardValues.SrCa = 8.838e-3;
PrimaryStandardValues.BaCa = 7.465e-6;
PrimaryStandardValues.UCa  = 1.192e-6;

% Global handle to UI Figure
UI = uifigure('Name','MATLAB Laser Reduction','WindowState','Maximized','Color',[0.3 0.3 0.3]);
drawnow; pause(1);

% Create tab groups:
tabdim = round([(UI.Position(3) - UI.Position(1))*0.1 (UI.Position(4) - UI.Position(2))*0.1 (UI.Position(3) - UI.Position(1))*0.8 (UI.Position(4) - UI.Position(2))*0.8]);
tabgp = uitabgroup(UI,'Position',tabdim,'TabLocation','Left');
tabnames = {'Load Data','Baseline Subtraction','Peak Identification','Measurement QA','Standard Correction','Align Laser Tracks','Match Data to Tracks'};
for i = 1:length(tabnames)
    tab{i} = uitab(tabgp,'Title',tabnames{i},'BackgroundColor','w','ForegroundColor','k');
end

    %% Part 1: Data Reduction
    
    % Load in all data in run:
    tabCount = UIprocesses(1);
    [RunDat,fnames] = Read_Data;
    
    % Identify and subtract baseline from the whole signal:
    tabCount = UIprocesses(tabCount);
    RunDatCor = Baseline_Subtraction(RunDat);
    
    % Prompt to see if peak IDs have already been made:
    tabCount = UIprocesses(tabCount);
    if sv == true
        % Load saved peaks from file:
        [file,path] = uigetfile([folder,filesep,'*.*'],'Open saved peak file');
        load(fullfile(path,file));
        figure(UI)
    else
        % Identify peaks by user input:
        delete(tabgp.SelectedTab.Children)
        [PeakRegions,PeakNames,PeakTimes,stdnames] = Peak_Identification(RunDatCor,fnames,tabCount);
    end
    
    % QA regions of peaks to be used, and average count ratios at spots:
    tabCount = UIprocesses(tabCount);
    [SpotData] = Peak_QA(PeakRegions,PeakNames,RunDatCor,stdnames);
    
    % Extract Element Ratios based on Standards:
    tabCount = UIprocesses(tabCount);
    SpotDatCor = StandardFit(SpotData,stdnames,PrimaryStandardValues);
    
    
    %% Part 2: Data mapping and plotting:

    % Load in image and define track paths:
    tabCount = UIprocesses(tabCount);
    [TrackCoords] = TrackAligner(tabCount-1);
    
    % Match Track coords with laser spot data:
    tabCount = UIprocesses(tabCount);
    [CombData] = DataCombiner(SpotDatCor,TrackCoords);
    
    