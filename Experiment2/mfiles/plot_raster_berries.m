function ber_raster = plot_raster_berries( data_subj, clrs)

hold all;

% berry_cell = cell(1, length(probe_indices));
idx_add = 0;
blk_add = 100;
for fn = 1:2
    
    berry_cell = [data_subj{fn}.berry_diff];
    probe_indices = data_subj{fn}.probe_indices;
    
    for trialCount =1:100+ blk_add % length(berry_cell) %length(probe_indices)
        
        trc = trialCount ;
%         if trc == 21 || trc == (31)
%            idx_add = idx_add + 5; 
% %            blk_add = trialCount + idx_add; 
%         end
%         
        
%         spikePos = [1:length(berry_cell{probe_indices(trialCount)})];%tVec(spikeMat(trialCount, :));
%         berry_arr = berry_cell{probe_indices(trialCount)};
        
        spikePos = [1:length(berry_cell{trialCount})];
        berry_arr = [0,berry_cell{trialCount}];% diff is one time index short
        berr_log = logical(berry_arr>0);
        spikePos = spikePos(berr_log);
        time_ = data_subj{fn}.time_ms{trialCount};
        patch_entry_ix = data_subj{fn}.enter_patch_ix(trialCount);
        patch_exit_ix = data_subj{fn}.leave_patch_ix(trialCount);
        trc
        for spikeCount = 1:length(spikePos)
            ber = berry_arr(spikePos(spikeCount));
            
            time_since_patch_entry = time_(spikePos(spikeCount)) -...
                time_(patch_entry_ix);
            plot([time_since_patch_entry time_since_patch_entry], ...
                idx_add + [trc-0.5 trc+0.5], ...
                clrs{fn},'LineWidth', ber*0.1);
            hold on
        end
        
        time_patch_exit = time_(patch_exit_ix) -...
                time_(patch_entry_ix);
        plot(time_patch_exit, idx_add + trc, 'k*','MarkerSize',5);
        hold on
    end
    %     idx_add = idx_add + trc + 5;
    xlabel('Time (s)');
    ylabel('Trial');
    % ylim([0 size(berry_cell, 1)+1]);
    
    beautifyfig
    keyboard
    figure;
    idx_add=0;
    
end
