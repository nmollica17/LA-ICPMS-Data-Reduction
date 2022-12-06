%% Laser Data processing

function [Outcoords] = TrackAligner(tabcount)
    global tab
    
    % Get colors for plotting:
    colors = get(groot,'DefaultAxesColorOrder');

    % Find and read in image:
    imfile = imgetfile;
    img = imread(imfile);

    % Display image on left:
    tabdim = tab{tabcount}.Position;
    ax1 = uiaxes(tab{tabcount},'BackgroundColor','w',...
        'Position',round([(tabdim(3) - tabdim(1))*0.1 (tabdim(4) - tabdim(2))*0.1...
                       (tabdim(3) - tabdim(1))*0.45 (tabdim(4) - tabdim(2))*0.8]),...
        'XTickLabel',[],'YTickLabel',[]);
    imagesc(ax1,img); axis(ax1,'image')
    
    [file,path] = uigetfile('*','Open exported tracks');
    
    % Load in coordinates of all points in a run, attach number ids to them
    [~,rawCoords] = xlsread([path,file]);
    coords = []; count = 1;
    for i = 1:length(rawCoords)
        if isempty(rawCoords{i,2})
            % read in and convert to mm:
            coords(count,:) = [count eval(rawCoords{i,8}).*1e-3];
            count = count + 1;
        end
    end
    
    % Plot coordinates on image, roughly line up with tracks
    ax2 = uiaxes(tab{tabcount},'BackgroundColor','w',...
    'Position',round([(tabdim(3) - tabdim(1))*0.55 (tabdim(4) - tabdim(2))*0.1...
                   (tabdim(3) - tabdim(1))*0.45 (tabdim(4) - tabdim(2))*0.8]),...
    'XTickLabel',[],'YTickLabel',[]);
    title(ax2,'Select groupings from TOP TO BOTTOM of coral')
    plot(ax2,coords(:,2),coords(:,3),'.k'); axis(ax2,'equal'); set(ax2,'ydir','reverse')
    text(ax2,coords(1,2),coords(1,3),'1','FontSize',20)
    xlim(ax2,[min(coords(:,2))*0.8 max(coords(:,2))*1.2])
    ylim(ax2,[min(coords(:,3))*0.8 max(coords(:,3))*1.2])
    hold(ax2,'on')
    
    % Prompt for user input for number of corals and number of groups:
    prompt = {'How many distinct corals are in the image?'};
    numcorals = inputdlg(prompt,'Corals',1); numcorals = str2double(numcorals{1});

    % Loop through corals and prompt user to ID points.
    for i = 1:numcorals
        count = 1;
        while numcorals >= 1
            roi = drawrectangle(ax2,'Color',colors(i,:));
            corners = roi.Vertices;
            
            % Select coords that are within the bounds:
            coordIDs = coords(:,2) > min(corners(:,1)) & coords(:,2) < max(corners(:,1))...
                     & coords(:,3) > min(corners(:,2)) & coords(:,3) < max(corners(:,2));
            
            % Write grouping to coords variable:
            coords(coordIDs,5) = i;
            coords(coordIDs,6) = count;
                 
            % Replot and delete rectangle:
            tempcrds = coords(coordIDs,:);
            plot(ax2,tempcrds(:,2),tempcrds(:,3),'.','color',colors(i,:))
            delete(roi)
            
            % Prompt to select region or continue:
            answer = questdlg(['Are there additional groupings?'],'Continue or move to next coral','More groupings','Next coral','Finished','Finished');
            if strcmp(answer,'Next coral')
                break;
            elseif strcmp(answer,'Finished')
                break;
            end  
            
            count = count +1;
        end
        
        % If ending early, remove unassigned points:
        if strcmp(answer,'Finished')
           coords(coords(:,6)==0,:) = [];
           break
        end
    end
    
    % Save figure with points added:
    names = strsplit(imfile,'.');
    imName = [names{1},'_IDs.','png'];
    % saveas(f,imName)
    
    %% Track Concatination
    
    % Define sections in grayscale for visual distinction:
    grayscale = repmat([0:100:250]',1,3)./255;
    
    % Initialize output array:
    Outcoords = cell(1,numcorals);
    
    % Concatinate each defined coral into one track (mm down core).
    for i = 1:numcorals
        numgroups = max(coords(coords(:,5)==i,6));
        % For each grouping, find axis and normalize axis to vertical:
        joinpt = 0;
        for j = 1:numgroups
            % Isolate coords for section:
            tempcrds = coords(coords(:,5)==i & coords(:,6)==j,:);
            % translate to depth down core:
            tempcrds(:,2:3) = tempcrds(:,2:3) - [min(tempcrds(:,2)) min(tempcrds(:,3))];
            
            % Fit major axis regression:
            [m,~] = lsqfitma(tempcrds(:,2),tempcrds(:,3));
            
            % Perform rotation:
            if abs(m) >= 1 % Slide is vertical
                theta = (90-atand(abs(m)))*sign(m);
                R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; %rotation matrix
                tempcrds(:,7:8) = [R*tempcrds(:,2:3)']';
                tempcrds(:,7:8) = tempcrds(:,7:8) - [mean(tempcrds(:,7)) min(tempcrds(:,8))] + joinpt;
                % save join point:
                joinpt = [mean(tempcrds(:,7)) max(tempcrds(:,8))];
            else % slide is horizontal
                theta = -atand(m)*sign(m);
                R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; %rotation matrix
                tempcrds(:,7:8) = [R*tempcrds(:,2:3)']';
                tempcrds(:,7:8) = tempcrds(:,7:8) - [min(tempcrds(:,7)) mean(tempcrds(:,8))];
                % determine up direction:
                if tempcrds(1,2) > tempcrds(end,2)
                    % if points are right to left, switch:
                    tempcrds(:,7) = abs(tempcrds(:,7)-max(tempcrds(:,7)));
                end
                % save join point:
                tempcrds(:,7:8) = tempcrds(:,7:8) + joinpt;
                joinpt = [max(tempcrds(:,7)) mean(tempcrds(:,8))];
            end
            

            
            % Plot points raw and adjusted:
            f = figure(1);
            subplot(121); hold off
            title(['Unrotated Coordinates - Coral ',num2str(i)]); hold on; axis equal;% set(gca,'YDir','Reverse')
            plot(coords(coords(:,5)==i & coords(:,6)==j,2),coords(coords(:,5)==i & coords(:,6)==j,3),'.','color',grayscale(j,:))
            subplot(122); hold off
            title(['Rotated Coordinates - Coral ',num2str(i)]); hold on; axis equal; set(gca,'YDir','Reverse')
            if abs(m) >= 1
                 plot(tempcrds(:,7),abs(tempcrds(:,8)),'.','color',grayscale(j,:))              
            else
                plot(tempcrds(:,8),abs(tempcrds(:,7)),'.','color',grayscale(j,:))
            end
            % save aligned coordinates to output,sorting by point ID:
            Outcoords{i} = sortrows(cat(1,Outcoords{i},[tempcrds(:,1) tempcrds(:,7) tempcrds(:,8)]),1);
        end
               
        % Pause, and save the figure:
        pause(2)
        saveas(f,[path,'Rotated Coords Coral ',num2str(i),'.png'])
        close
    end   
    
end
