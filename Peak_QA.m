%% Laser Data processing

function [SpotDat] = Peak_QA(PeakRegions,PeakNames,RawData,stdnames)
    
    global UI
    
    d = uiprogressdlg(UI,'Title','Data Quality Control',...
        'Message','Identifying suitable peak regions for analysis');
    
    fields = fieldnames(RawData);
    
    % Pull out and evaluate signal related to each peak ID
    countNum = 1;
    for i = 1:length(RawData.time)
        count1 = 1; count2 = 1;
        peakNums = 1:max(PeakRegions{i});
        for j = 1:max(peakNums)
            % Create seperate temp variable to keep track of time:
            t  = RawData.time{i}(PeakRegions{i}==j); 
            % Must have a minimum of 20 seconds of sample:
            if max(t) - min(t) > 20
                [~,eid] = min(abs(t-t(1)-50)); % closest value to 50 (end buffer).
                
                % clip data to cutoff window
                for k = 1:length(fields)
                    peakData.(fields{k}) = RawData.(fields{k}){i}(PeakRegions{i}==j);
                    peakData.(fields{k}) = peakData.(fields{k})(t>=t(1)+5 & t<=t(eid)-5);
                end
                                
                % Second, perform multivariate cluster analysis using Sr/Ca/U
                nsmooth = 5; % number of points smoothed in cluster analysis
                nclusters = 3; % number of possible major clusters (?)
                clusters = clusterdata([smooth(peakData.Sr88,nsmooth) smooth(peakData.U238,nsmooth) smooth(peakData.Ca43,nsmooth)],'Maxclust',nclusters,'distance','seuclidean');

                % Identify minor clusters:

                isminor = false(size(clusters));
                for c = 1:max(clusters)
                    if round( (sum(clusters==c)/length(clusters))*nclusters ) ==0
                        isminor(clusters==c) = true;
                    end
                end
                
                % Remove minor clusters from data to be analyzed:
                clusters(isminor) = NaN;
                [largestCluster,lcLength] = mode(clusters(~isnan(clusters)));
                
                % test if means are the same if there are 2 or more clusters:
                if length(unique(clusters(~isnan(clusters)))) > 1
                    [~,~,SrStats] = anova1(peakData.Sr88./peakData.Ca43,clusters,'off'); SrStats = multcompare(SrStats,'Display','off');
                    [~,~,UStats]  = anova1(peakData.U238./peakData.Ca43,clusters,'off');  UStats = multcompare(UStats,'Display','off');
                    
                    % Pick out pairs that are NOT different in mean Sr/Ca and U/Ca value:
                    successes = SrStats(SrStats(:,6) > 0.05 & UStats(:,6) > 0.05,1:2);
                    
                    % See if any pairing is larger than the largest cluster, otherwise use the largest cluster
                    
                    for k = 1:size(successes,1)
                        testPair = clusters(clusters == successes(k,1) | clusters == successes(k,2));
                        if length(testPair) > lcLength
                           largestCluster = successes(k,:);
                           lcLength = length(testPair); % update length of largest cluster
                        end
                    end
                end
                
                % kill all but the largest cluster (or pair)
                for k = 1:length(fields)
                    if length(largestCluster) == 1
                        peakData.(fields{k}) = peakData.(fields{k})(clusters == largestCluster);
                    else
                        peakData.(fields{k}) = peakData.(fields{k})(clusters == largestCluster(1) | clusters == largestCluster(2));
                    end
                end

                
                % Convert to element ratios:
                enames = fields(contains(fields,digitsPattern)&~contains(fields,'Ca43'));
                
                % Split points into primary, secondary, sample, other:
                if strcmp(PeakNames{i}{j},stdnames.Primary)
                    SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).t(count1) = mean(peakData.time);
                    SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).n(count1) = numel(peakData.time);
                    for k = 1:length(enames)
                        SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).([enames{k}(isstrprop(enames{k},'alpha')),'Ca'])(count1) ...
                                 = median(peakData.(enames{k})./peakData.Ca43);
                        SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).([enames{k}(isstrprop(enames{k},'alpha')),'Cae'])(count1) ...
                                 = std(peakData.(enames{k})./peakData.Ca43);
                    end
                    count1 = count1 +1;
                elseif strcmp(PeakNames{i}{j},stdnames.Secondary)
                    SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).t(count2) = mean(peakData.time);
                    SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).n(count2) = numel(peakData.time);
                    for k = 1:length(enames)
                        SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).([enames{k}(isstrprop(enames{k},'alpha')),'Ca'])(count2) ...
                                 = median(peakData.(enames{k})./peakData.Ca43);
                        SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).([enames{k}(isstrprop(enames{k},'alpha')),'Cae'])(count2) ...
                                 = std(peakData.(enames{k})./peakData.Ca43);
                    end
                    count2 = count2 +1;
                else
                    try
                        SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).t = mean(peakData.time);
                        SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).n = numel(peakData.time);
                        for k = 1:length(enames)
                            SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).([enames{k}(isstrprop(enames{k},'alpha')),'Ca']) ...
                                     = median(peakData.(enames{k})./peakData.Ca43);
                            SpotDat.(['Run',num2str(i)]).(PeakNames{i}{j}).([enames{k}(isstrprop(enames{k},'alpha')),'Cae']) ...
                                     = std(peakData.(enames{k})./peakData.Ca43);
                        end
                    catch
                        % Add letter in front of numbers:
                        SpotDat.(['Run',num2str(i)]).(['S_',num2str(countNum)]).t = mean(peakData.time);
                        SpotDat.(['Run',num2str(i)]).(['S_',num2str(countNum)]).n  = numel(peakData.time);
                        for k = 1:length(enames)
                            SpotDat.(['Run',num2str(i)]).(['S_',num2str(countNum)]).([enames{k}(isstrprop(enames{k},'alpha')),'Ca']) ...
                                     = median(peakData.(enames{k})./peakData.Ca43);
                            SpotDat.(['Run',num2str(i)]).(['S_',num2str(countNum)]).([enames{k}(isstrprop(enames{k},'alpha')),'Cae']) ...
                                     = std(peakData.(enames{k})./peakData.Ca43);
                        end
                        countNum = countNum + 1;
                    end
                end

            else
                % Still increment count Num for skipped sample peaks:
                if ~isnan(str2double(PeakNames{i}{j}))
                    countNum = countNum + 1;
                end
            end
        end
        
        d.Value = min(d.Value + 1/length(RawData.time),1);
    end
end
