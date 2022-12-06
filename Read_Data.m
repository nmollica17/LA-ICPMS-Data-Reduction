%% Laser Data processing

function [RunData,fnames] = Read_Data
    global UI folder
    
    % Prompt user to select laser data sheets:
    [fnames,folder] = uigetfile('*','Open laser data sheets','multiselect','on');
    fnames = fullfile(folder,fnames)';
    figure(UI)
    
    if ischar(fnames)
        fnames = cellstr(fnames);
    end
    
    % setup progess dialog:
    d = uiprogressdlg(UI,'Title','Reading Data',...
        'Message','Please wait while the selected files are loaded.');
    
    % Load in raw data:
    for i = length(fnames):-1:1
        rawdat = xlsread(fnames{i});
        RunData.time{i} = rawdat(1:end,1);
        RunData.B11{i}  = rawdat(1:end,2);
        RunData.Mg25{i} = rawdat(1:end,3);
        RunData.Ca43{i} = rawdat(1:end,4);
        RunData.Sr88{i} = rawdat(1:end,5);
        RunData.Ba138{i}= rawdat(1:end,6);
        RunData.U238{i} = rawdat(1:end,7);
        RunData.TotalBeam{i} = sum(rawdat(1:end,2:7),2);
        d.Value = min(d.Value + 1/length(fnames),1);
    end
end