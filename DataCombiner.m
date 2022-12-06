%% Combine structured output into large matrix:

function [CombData] = DataCombiner(SpotDatCor,TrackCoords)

    % Loop across corals:
    for c = 1:length(TrackCoords)
       
        % Read out run names:
        runs = fieldnames(SpotDatCor);
        
        spotnames = fieldnames(SpotDatCor.(runs{1}));
        
        % Loop down track coordinates and index out SpotData
        runCounter = 1; j = 1; spotnames = fieldnames(SpotDatCor.(runs{runCounter}));
        while j <= length(TrackCoords{c})
            index = TrackCoords{c}(j,1);
            
            %loop down spotnames to the end of the run data:
            for k = 1:length(spotnames)
                if index == str2double(erase(spotnames{k},'S_'))
                    % Concatinate data, taking absolute value of coordinates
                    CombData{c}(j,:) = [TrackCoords{c}(j,1:3),...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).SrCa_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).SrCae_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).UCa_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).UCae_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).BCa_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).BCae_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).BaCa_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).BaCae_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).MgCa_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).MgCae_m, ...
                        SpotDatCor.(runs{runCounter}).(spotnames{k}).n];
                        j = j + 1;
                    break
                elseif index < str2double(erase(spotnames{k},'S_'))
                    % The point was removed, and we have passed it
                    j = j + 1;
                    break
                end
            end
            
            % If the spot is not found, or is the last in the run, move to next run
            if k==length(spotnames) && runCounter < length(runs)
            	runCounter = runCounter + 1;
                spotnames = fieldnames(SpotDatCor.(runs{runCounter}));
            elseif k == length(spotnames) && runCounter == length(runs)
                break
            end
        end
    end
end
