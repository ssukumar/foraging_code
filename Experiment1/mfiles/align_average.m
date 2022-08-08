function avg_vel =  align_average( vel_profile_cell,max_ix_mat,c,...
    probe_nonprobe, probe_indices)
% Gets a cell of time series velocity profiles
% Aligns them to the peak and then  averages
%
%probe_nonprobe: flag to determine whether to separate probe or non probe


num_trials = length(vel_profile_cell);
vel_profile_mat = [];

for tr = 2:num_trials
    tr1 = tr-1;
    ix_max = max_ix_mat(tr);
    
    if isnan(ix_max) && tr1==1
        keyboard
    end
    %     tr
    if ~isnan(ix_max) && ~isempty(vel_profile_cell{tr})
        
        if vel_profile_cell{tr}(max_ix_mat(tr)) < 0
            vel = -1*(vel_profile_cell{tr});
        else
            vel = vel_profile_cell{tr};
            
        end
        
        if tr >2
            
            if ix_max > prev_ix_max
                
                temp= [nan(tr1-1,(ix_max - prev_ix_max)),...
                    vel_profile_mat(1:tr1-1, :)];
                vel_profile_mat = temp;
                prev_ix_max = ix_max;
                
                if size(vel_profile_mat,2) > length(vel)
                    vel_fill = [vel, nan(1, size(temp,2)-length(vel))];
                else
                    vel_fill = vel;
                    temp = [vel_profile_mat(1:tr1-1,:), nan(tr1-1, length(vel)-size(temp,2))];
                    vel_profile_mat = temp;
                end
                vel_profile_mat(tr1,:) = vel_fill;
                
            else
                
                vel_fill = [nan(1,prev_ix_max - ix_max) , vel];
                if length(vel_fill) > size(vel_profile_mat,2)
                    
                    temp = [vel_profile_mat(1:tr1-1, :), nan(tr1-1, length(vel_fill)- ...
                        size(vel_profile_mat,2))];
                    vel_profile_mat = temp;
                    vel_profile_mat(tr1,:) = vel_fill;
                    
                else
                    vel_fill = [vel_fill, nan(1,size(vel_profile_mat,2) -...
                        length(vel_fill) )];
                    vel_profile_mat(tr1, :) = vel_fill;
                end
                
            end
            
        else
            vel_profile_mat(tr1,:) = vel_profile_cell{tr};
            prev_ix_max = ix_max;
        end
        
    else
        vel_profile_mat(tr1,:) = nan;
        
    end
    %     prev_ix_max = ix_max;
    %     vel_profile_mat
    %     keyboard
    
    
    
end
 
% 
% for tr1=1:num_trials-1
%     plot(vel_profile_mat(tr1,:))
%     drawnow
%     keyboard
%     hold on
% %     debugging only
%     title(sprintf('Trial # %d',tr1))
% 
% end
% 
% keyboard

if probe_nonprobe==0
    
    avg_vel_prof = nanmean(vel_profile_mat(2:end,:),1);
    ster_vel_prof = nanstd(vel_profile_mat(2:end,:),0,1);
    
    idx = find(ster_vel_prof,1,'first');
    avg_vel_prof = avg_vel_prof(idx:end);
    ster_vel_prof = ster_vel_prof(idx:end);
    keepIndex = ~isnan(avg_vel_prof) & ~isnan(ster_vel_prof);
    avg_vel_prof = avg_vel_prof(keepIndex);
    ster_vel_prof = ster_vel_prof(keepIndex);
    
    gcf;
    
    if ~isempty(avg_vel_prof)
        plotsdshading([1:length(avg_vel_prof)],avg_vel_prof',ster_vel_prof', c)
        xlabel('Time (s)','Fontsize',13)
        ylabel('Velocity (m/s)','Fontsize',13);
    end
    
    avg_vel = avg_vel_prof;
    
else
    avg_vel = cell(1,2);
    p_indices = {setxor(probe_indices,2:size(vel_profile_mat,1)),...
        probe_indices};
    for p = 1:2
        indices = p_indices{p};
        
        if p==1
            clr = c;
        else
            clr='g';
        end
        
        indices = indices(indices<size(vel_profile_mat,1));
        
        avg_vel_prof = nanmean(vel_profile_mat(indices,:),1);
        ster_vel_prof = nanstd(vel_profile_mat(indices,:),0,1)./...
            sqrt(length(vel_profile_mat(indices,:))); % standard error of mean
        
        idx = find(ster_vel_prof,1,'first');
        avg_vel_prof = avg_vel_prof(idx:end);
        ster_vel_prof = ster_vel_prof(idx:end);
        keepIndex = ~isnan(avg_vel_prof) & ~isnan(ster_vel_prof);
        avg_vel_prof = avg_vel_prof(keepIndex);
        ster_vel_prof = ster_vel_prof(keepIndex);
        
        gcf;
        
        if ~isempty(avg_vel_prof)
            plotsdshading([1:length(avg_vel_prof)],avg_vel_prof',...
                ster_vel_prof', clr)
%             xlabel('Time (s)','Fontsize',13)
%             ylabel('Grip Force (N)','Fontsize',13);
        end
        
        %         end
        avg_vel{p} = avg_vel_prof; 
    end
end

