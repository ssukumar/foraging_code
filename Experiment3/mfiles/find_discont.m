%% Find Discontinuities and fix them.
% Updated 5-10-2019 by Gary Bruening
% This function takes in the x-y position data and uses the average changes
% in movement to see if there are spikes in the data that needs to be
% removed.
%
% Once the errors are found, use spline interpolation to fill in the spiked
% data.

% This function is hugely important for peak velocity calculation as the
% spikes can be REAL big.

function Data = find_discont(Px, Py, V, Data, c,subj)

% assumes each trial is a col
% V can be Vy or Vx

scrsz = get(0,'ScreenSize');

% figure, plot(range(diff(V)),'o'), ylabel('range of diff(V)'), xlabel('trial')
Data.discont.trials = find(range(diff(V)) > mean(range(diff(V)))*4);

for i = 1:length(Data.discont.trials)
    
    Data.discont.nspikes(i)=length(find(diff(V(:,Data.discont.trials(i)))>mean(range(diff(V)))*2));
    
    for k=1:Data.discont.nspikes(i)
        if isempty(find(diff(Data.vy(:,Data.discont.trials(i)))>mean(range(diff(V)))*2))
            break
        end

        Data.discont.frames = find(diff(V(:,Data.discont.trials(i)))>mean(range(diff(V)))*2,k);
        Data.discont.refframe(i,1) = Data.discont.frames(k);
            
        Data.discont.frame(i,1)=Data.discont.refframe(i,1)-1;
        Data.discont.frame(i,2)=Data.discont.refframe(i,1)+1;

        % Sometimes the spike happens at the beginning/end of trial and it
        % can't find other side, so it picks a frame at the max of the movement
        % and results in beginning/end frames that are out of order.
        if Data.discont.frame(i,1)>Data.discont.frame(i,2)
            num_frames = sum(~isnan(V(:,Data.discont.trials(i))));
            if abs(Data.discont.frame(i,1) - num_frames)<=3
                Data.discont.frame(i,2) = Data.discont.frame(i,1); % at end of mvt
            elseif Data.discont.frame(i,2) <=3
                Data.discont.frame(i,1) = Data.discont.frame(i,2); % at start of mvt
            else
                fprintf('\n Frames are out of order: trial %i \n',i)
            end
        end

            dt = Data.discont.trials(i);
            df = round(mean(Data.discont.frame(i,1:2)));
            Data.x(df,dt) = mean(Data.x(Data.discont.frame(i,:),dt));
            Data.y(df,dt) = mean(Data.y(Data.discont.frame(i,:),dt));

        vlength=sum(~isnan(Data.vy(:,Data.discont.trials(i))));    

        if (Data.discont.frame(i)-vlength) < -3
            xi = Data.discont.frame(i,1)-2:Data.discont.frame(i,2)+2;
            Data.x(xi,dt)=nan;Data.y(xi,dt)=nan;
            Data.x(xi,dt)=spline(1:length(Data.x(:,dt)),Data.x(:,dt),xi);
            [a,b]=lastwarn;
            warning('off',b);
            Data.y(xi,dt)=spline(1:length(Data.y(:,dt)),Data.y(:,dt),xi);

            Data.vx(xi,dt)=nan;Data.vy(xi,dt)=nan;
            Data.vx(xi,dt)=spline(1:length(Data.vx(:,dt)),Data.vx(:,dt),xi);                
            Data.vy(xi,dt)=spline(1:length(Data.vy(:,dt)),Data.vy(:,dt),xi);

            Data.fx(xi,dt)=nan;Data.fy(xi,dt)=nan;
            Data.fx(xi,dt)=spline(1:length(Data.fx(:,dt)),Data.fx(:,dt),xi);
            Data.fy(xi,dt)=spline(1:length(Data.fy(:,dt)),Data.fy(:,dt),xi);

        elseif (Data.discont.frame(i)-vlength) <= -1
            xi = Data.discont.frame(i,1)-1:Data.discont.frame(i,2)+1;
            
            Data.x(xi,dt) = interp1([xi(1) xi(end)], Data.x([xi(1) xi(end)],dt), xi);                
            Data.y(xi,dt) = interp1([xi(1) xi(end)], Data.y([xi(1) xi(end)],dt), xi);

            Data.vx(xi,dt) = interp1([xi(1) xi(end)], Data.vx([xi(1) xi(end)],dt), xi);
            Data.vy(xi,dt) = interp1([xi(1) xi(end)], Data.vy([xi(1) xi(end)],dt), xi);

            Data.fx(xi,dt) = interp1([xi(1) xi(end)], Data.fx([xi(1) xi(end)],dt), xi);
            Data.fy(xi,dt) = interp1([xi(1) xi(end)], Data.fy([xi(1) xi(end)],dt), xi);

        end            
    end
end

