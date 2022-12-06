%% Laser Data processing

function [SpotData] = StandardFit(SpotData,stdnames,psv)

    global UI
    
    d = uiprogressdlg(UI,'Title','Standard Correction',...
        'Message','Converting to molar ratios');

    rfields = fieldnames(SpotData);
    ERnames = fieldnames(SpotData.(rfields{1}).(stdnames.Primary));
    ERnames = ERnames(contains(ERnames,'Ca')&~contains(ERnames,'Cae'));
    
    % Iterate through runs and fit various functions to the standards
    for i = 1:length(rfields)
        % Unpack standard points:
        t = SpotData.(rfields{i}).(stdnames.Primary).t;

        % Unpack Query (sample) points:
        sfields = fieldnames(SpotData.(rfields{i})); count = 1;
        for j = 1:length(sfields)
            % Only do samples
            if contains(sfields{j},'S_')
                tq(count) = SpotData.(rfields{i}).(sfields{j}).t;
                count = count + 1;
            end
        end
        
        % Bootstrap interpolation using algorithm of choice (here I use pchip):
        n = 1e5; Primary_fit = zeros(length(tq),length(ERnames),n);
        for k = 1:length(ERnames)
            parfor j = 1:n
                % Simulate data
                tempER = normrnd(SpotData.(rfields{i}).(stdnames.Primary).(ERnames{k}),...
                                 SpotData.(rfields{i}).(stdnames.Primary).([(ERnames{k}),'e']));
                % Calculate fit
                Primary_fit(:,k,j) = pchip(t,tempER,tq);
            end
        end
        
        %Average simulated fits:
        Primary_fit = mean(Primary_fit,3);
        
        % Convert to molar ratios:
        count = 1;
        for j = 1:length(sfields)
            % Only do samples
            if contains(sfields{j},'S_')
                for k = 1:length(ERnames)
                    SpotData.(rfields{i}).(sfields{j}).([(ERnames{k}),'_m']) ...
                        = SpotData.(rfields{i}).(sfields{j}).(ERnames{k})/Primary_fit(count,k)*psv.(ERnames{k});
                    SpotData.(rfields{i}).(sfields{j}).([(ERnames{k}),'e_m']) ...
                        = SpotData.(rfields{i}).(sfields{j}).([(ERnames{k}),'e'])/Primary_fit(count,k)*psv.(ERnames{k});
                end
                count = count + 1;
            end
        end
        
        % UNDER REVISIONS
%         % unpack and convert secondary standard:
%         t2 = SpotData.(rfields{i}).(stdnames.Secondary).t;
%         SrCa2 = SpotData.(rfields{i}).(stdnames.Secondary).SrCa;
%         UCa2  = SpotData.(rfields{i}).(stdnames.Secondary).UCa;
%         
%         % Fit using same algorithm:
%         SrCa_fit2 = pchip(t,SrCa,t2);
%         UCa_fit2  = pchip(t,UCa,t2);
%         
%         %Save back into structure:
%         SpotData.(rfields{i}).(stdnames.Secondary).SrCa_m = SpotData.(rfields{i}).(stdnames.Secondary).SrCa./SrCa_fit2.*psv.SrCa;
%         SpotData.(rfields{i}).(stdnames.Secondary).UCa_m = SpotData.(rfields{i}).(stdnames.Secondary).UCa./UCa_fit2.*psv.UCa;
    
        d.Value = min(d.Value + 1/length(rfields),1);
    end  
end
