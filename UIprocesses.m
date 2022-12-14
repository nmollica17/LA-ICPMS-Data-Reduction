%% Notification Boxes for Laser Processor

function count = UIprocesses(count)
    global UI tabgp tab
    
    tabdim = tab{1}.Position;
    
    if count == 1 % Load Data
        % In this tab, display the welcome message and have a load data button
        uilabel(tab{count},'Text',sprintf('First, you will need to load the raw laser datafiles.'),...
                     'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                                       (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
                      'FontSize',18,...
                      'HorizontalAlignment','center');
        uibutton(tab{count},'Text','Load Some Data','ButtonPushedFcn','uiresume(UI)',...
                     'Position',round([(tabdim(3) - tabdim(1))*0.35 (tabdim(4) - tabdim(2))*0.7...
                                       (tabdim(3) - tabdim(1))*0.15 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{count},'Text','Quit MATLAB','ButtonPushedFcn','uiresume(UI);quit',...
                     'Position',round([(tabdim(3) - tabdim(1))*0.55 (tabdim(4) - tabdim(2))*0.7 ...
                                       (tabdim(3) - tabdim(1))*0.15 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uiwait(UI)
        count = 2;

    elseif count == 2 % Automated Background Subtraction
        tabgp.SelectedTab = tab{count}; drawnow;
        uilabel(tab{count},'Text',sprintf('Background removal does not have options right now, sorry!'),...
             'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                               (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
              'FontSize',18,...
              'HorizontalAlignment','center');
        count = 3; 
        
    elseif count == 3 % Peak Identification tab
        tabgp.SelectedTab = tab{count}; drawnow;
        uilabel(tab{count},'Text',sprintf('Have you already identified peaks for this data?'),...
             'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                               (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
              'FontSize',18,...
              'HorizontalAlignment','center');
        uibutton(tab{count},'Text','I have Saved Peaks','ButtonPushedFcn','sv = true;uiresume(UI)',...
                     'Position',round([(tabdim(3) - tabdim(1))*0.35 (tabdim(4) - tabdim(2))*0.7...
                                       (tabdim(3) - tabdim(1))*0.15 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uibutton(tab{count},'Text','I will ID new Peaks','ButtonPushedFcn','sv = false;uiresume(UI)',...
                     'Position',round([(tabdim(3) - tabdim(1))*0.55 (tabdim(4) - tabdim(2))*0.7 ...
                                       (tabdim(3) - tabdim(1))*0.15 (tabdim(4) - tabdim(2))*0.05]),...
                      'FontSize',14,...
                      'HorizontalAlignment','center');
        uiwait(UI)
        count = 4;
        
    elseif count == 4 % Data QA
        tabgp.SelectedTab = tab{count}; drawnow;
        uilabel(tab{count},'Text',sprintf('Data QA also does not have options right now, sorry!'),...
             'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                               (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
              'FontSize',18,...
              'HorizontalAlignment','center');
        count = 5;
        
    elseif count == 5 % Standard Correction
        tabgp.SelectedTab = tab{count}; drawnow;
        uilabel(tab{count},'Text',sprintf('Standard Correction Options'),...
             'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                               (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
              'FontSize',18,...
              'HorizontalAlignment','center');
        count = 6;
        
    elseif count == 6 % Align laser tracks
        tabgp.SelectedTab = tab{count}; drawnow;
        uilabel(tab{count},'Text',sprintf('Track Alignment'),...
             'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                               (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
              'FontSize',18,...
              'HorizontalAlignment','center');
        count = 7;
        
    elseif count == 7 % Match data to tracks
        tabgp.SelectedTab = tab{count}; drawnow;
        uilabel(tab{count},'Text',sprintf('Track Assignment'),...
             'Position',round([(tabdim(3) - tabdim(1))*0.2 (tabdim(4) - tabdim(2))*0.8...
                               (tabdim(3) - tabdim(1))*0.6 (tabdim(4) - tabdim(2))*0.1]),...
              'FontSize',18,...
              'HorizontalAlignment','center');
        count = 8;               
    end
end