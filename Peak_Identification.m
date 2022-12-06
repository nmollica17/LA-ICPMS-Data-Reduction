%% Laser Data processing

function [PeakRegions,PeakNames,PeakTimes,stdnames] = Peak_Identification(RunData,fnames,tabcount)
    global UI tab t x PPMedianTime PeakPeriods PPNumbers bs stdInfo
    fields = fieldnames(RunData);
    i = 1;

    % Using background corrected data, identify peaks
    while i <= length(RunData.time)
        
        % Identify peaks by Ca43 counts
        PeakPeriods = RunData.Ca43{i} > 3e4;
        
        % Discard areas less than 10 seconds long
        PPThreshTime = 10;
        tstep = 0.5; % Sampling frequency is half a second on our machine
        PeakPeriods = bwareaopen(PeakPeriods,(PPThreshTime/tstep));
        
        % number each peak
        PPNumbers = bwlabel(PeakPeriods);
        stats = regionprops(PeakPeriods);
        PPMedianTime = RunData.time{i}(round(cat(1,stats.Centroid)));      
        
        % Prompt for standard spacing:
        prompt = {'Primary Standard Name','Spacing','First position in Run',...
                  'Secondary Standard Name','Spacing','First position in Run',...
                  'Bracketing Standards (list all separeted by spaces, i.e. "NIST NIST OTO..."'};
        num_lines = 1;
        def = {'JCP','11','8','JAR','11','7','NIST NIST OTO OTO BC BC'};
        stdInfo = inputdlg(prompt,'Standard Spacing',num_lines,def);
        bs = strsplit(stdInfo{7});
        stdnames.Primary = stdInfo{1}; stdnames.Secondary = stdInfo{4}; stdnames.Other = stdInfo{7};
        
        % Run Peak Labeler:
        [PeakLabels] = PeakLabeler(PPNumbers,bs,stdInfo);
        
        % Plot data, prompt for user input
        tabdim = tab{tabcount-1}.Position;
        ax = uiaxes(tab{tabcount-1},'BackgroundColor','w',...
                'Position',round([(tabdim(3) - tabdim(1))*0.1 (tabdim(4) - tabdim(2))*0.1...
                               (tabdim(3) - tabdim(1))*0.9 (tabdim(4) - tabdim(2))*0.8]),...
                'FontSize',18);
        ax.Interactions = [panInteraction zoomInteraction];
        plot(ax,RunData.time{i},RunData.TotalBeam{i});
        t = text(ax,PPMedianTime(:,2),zeros(size(PPMedianTime(:,2))),PeakLabels,'HorizontalAlignment','center','clipping','on',...
                 'ButtonDownFcn',@SaveXY);
        title(ax,'Click labels to change, or click Add..')
        xlimits = [0 1.5e3]; xlim(ax,xlimits); ylim(ax,[-0.5e7 max(RunData.TotalBeam{i})])
        
        % Place buttons for assignment:
        btn = uibutton(tab{3},'Text','Remove','ButtonPushedFcn',@(btn,event) RemoveFcn(RunData,ax,i),...
                     'Position',round([(tabdim(3) - tabdim(1))*1 (tabdim(4) - tabdim(2))*0.7...
                                       (tabdim(3) - tabdim(1))*0.13 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{3},'Text','Rename','ButtonPushedFcn',@(btn,event) RenameFcn,...
                     'Position',round([(tabdim(3) - tabdim(1))*1 (tabdim(4) - tabdim(2))*0.6...
                                       (tabdim(3) - tabdim(1))*0.13 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{3},'Text','Add...','ButtonPushedFcn',@(btn,event) AddFcn(RunData,ax,i),...
                     'Position',round([(tabdim(3) - tabdim(1))*1 (tabdim(4) - tabdim(2))*0.5...
                                       (tabdim(3) - tabdim(1))*0.13 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{3},'Text','Finished','ButtonPushedFcn',@(btn,event) delete(ax),...
                     'Position',round([(tabdim(3) - tabdim(1))*1 (tabdim(4) - tabdim(2))*0.4...
                                       (tabdim(3) - tabdim(1))*0.13 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{3},'Text',char(8592),'ButtonPushedFcn',@(btn,event) backFcn(ax),...
                     'Position',round([(tabdim(3) - tabdim(1))*1 (tabdim(4) - tabdim(2))*0.3...
                                       (tabdim(3) - tabdim(1))*0.06 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{3},'Text',char(8594),'ButtonPushedFcn',@(btn,event) forewardFcn(ax,RunData,i),...
                     'Position',round([(tabdim(3) - tabdim(1))*1.07 (tabdim(4) - tabdim(2))*0.3...
                                       (tabdim(3) - tabdim(1))*0.06 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{3},'Text','Redo Assgn','ButtonPushedFcn',@(btn,event) RedoFcn(ax),...
                     'Position',round([(tabdim(3) - tabdim(1))*1 (tabdim(4) - tabdim(2))*0.2...
                                       (tabdim(3) - tabdim(1))*0.13 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        waitfor(ax)
        [PeakLabels] = PeakLabeler(PPNumbers,bs,stdInfo);
        PeakRegions{i} = PPNumbers;
        PeakNames{i} = PeakLabels;
        PeakTimes{i} = PPMedianTime;
        i = i+1;
    end
    
    % Save peaks into .mat file
    [file,path] = uiputfile('*.mat','Save peak IDs');
    save(fullfile(path,file),'PeakRegions','PeakNames','PeakTimes','stdnames')  
end

    function SaveXY(src,~)
        global x
        if src.Color == [0 0 0]
            x = src.Position(1);
            src.Color = [1 0 0];
            src.FontWeight = 'Bold';
        else
            src.Color = [0 0 0];
            x = [];
            src.FontWeight = 'Normal';
        end
    end        

    function RemoveFcn(RunData,ax,i)
        global t x PPMedianTime PeakPeriods PPNumbers bs stdInfo
        % find and remove identified peak:
        [~,xid] = min(abs(PPMedianTime(:,2)-x));
        PeakPeriods(PPNumbers == xid) = 0;
        PPNumbers = bwlabel(PeakPeriods);
        stats = regionprops(PeakPeriods);
        PPMedianTime = RunData.time{i}(round(cat(1,stats.Centroid)));

        % Rerun peak labeling:
        [PeakLabels] = PeakLabeler(PPNumbers,bs,stdInfo);

        % Update plot labels
        delete(t);
        t = text(ax,PPMedianTime(:,2),zeros(size(PPMedianTime(:,2))),PeakLabels,'HorizontalAlignment','center','clipping','on',...
                 'ButtonDownFcn',@SaveXY);
    end

    function RenameFcn
        % Does nothing for right now
    end

    function AddFcn(RunData,ax,i)
        global t PPMedianTime PeakPeriods PPNumbers bs stdInfo
        
        % Give instructions:
        Options.Interpreter = 'tex';
        Options.WindowStyle = 'modal';
        uiwait(msgbox('\fontsize{10} Select peak area','Add new peak',Options))
        
        % Query new point from graph:
        roi = drawrectangle(ax);
        x = unique(roi.Vertices(:,1));
        delete(roi)
        
        % Insert new peak:
        [~,xid(1)] = min(abs(RunData.time{i}-x(1)));
        [~,xid(2)] = min(abs(RunData.time{i}-x(2)));
        PeakPeriods(xid(1):xid(2)) = 1;
        PPNumbers = bwlabel(PeakPeriods);
        stats = regionprops(PeakPeriods);
        PPMedianTime = RunData.time{i}(round(cat(1,stats.Centroid)));

        % Rerun peak labeling:
        [PeakLabels] = PeakLabeler(PPNumbers,bs,stdInfo);

        % Update plot labels
        delete(t);
        t = text(ax,PPMedianTime(:,2),zeros(size(PPMedianTime(:,2))),PeakLabels,'HorizontalAlignment','center','clipping','on',...
                 'ButtonDownFcn',@SaveXY);
    end  

    function backFcn(ax)
        currentX = xlim(ax);
        xlim(ax,currentX - min(1.5e3,currentX(1)));
    end
    
    function forewardFcn(ax,RunData,i)
        currentX = xlim(ax);
        xlim(ax,currentX + min(1.5e3,max(RunData.time{i}) - currentX(2)));
    end
    
    function RedoFcn(ax)
        global PPNumbers PPMedianTime t bs stdInfo
        % Reprompt for spacing info:
        prompt = {'Primary Standard Name','Spacing','First position in Run',...
          'Secondary Standard Name','Spacing','First position in Run',...
          'Bracketing Standards (list all separeted by spaces, i.e. "NIST NIST OTO..."'};
        num_lines = 1;
        def = {'JCP','11','8','JAR','11','7','NIST NIST OTO OTO BC BC'};
        stdInfo = inputdlg(prompt,'Standard Spacing',num_lines,def);
        bs = strsplit(stdInfo{7});
       
        % Relabel peaks:
        [PeakLabels] = PeakLabeler(PPNumbers,bs,stdInfo);
        
        % Update plot labels
        delete(t);
        t = text(ax,PPMedianTime(:,2),zeros(size(PPMedianTime(:,2))),PeakLabels,'HorizontalAlignment','center','clipping','on',...
                 'ButtonDownFcn',@SaveXY);
    end
    
function [PeakLabels] = PeakLabeler(PPNumbers,bs,stdInfo)
        % Create Peak Labels
        count = 1; PeakLabels = cell(max(PPNumbers),1);
        for k = 1:max(PPNumbers)
            if k <= length(strsplit(stdInfo{7})) && ~isempty(bs{1})
                PeakLabels{k} = bs{k};
            elseif k > max(PPNumbers) - length(strsplit(stdInfo{7})) && ~isempty(bs{1})
                try
                    PeakLabels{k} = bs{k - (max(PPNumbers)- length(strsplit(stdInfo{7})))};
                catch; break; end
            elseif k == str2double(stdInfo{3}) || mod(k-str2double(stdInfo{3})+1,str2double(stdInfo{2})+1) == 1
                PeakLabels{k} = stdInfo{1};
            elseif k == str2double(stdInfo{6}) || mod(k-str2double(stdInfo{6})+1,str2double(stdInfo{5})+1) == 1
                PeakLabels{k} = stdInfo{4};
            else
                PeakLabels{k} = [num2str(count)];
                count = count + 1;
            end
        end
end

