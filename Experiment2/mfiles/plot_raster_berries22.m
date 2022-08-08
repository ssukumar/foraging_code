function ber_raster = plot_raster_berries22( berry_cell, probe_indices, c)

hold all;

% berry_cell = cell(1, length(probe_indices));

for trialCount = 1:length(probe_indices)
    spikePos = [1:length(berry_cell{probe_indices(trialCount)})];%tVec(spikeMat(trialCount, :));
    berry_arr = berry_cell{probe_indices(trialCount)};
    berr_log = logical(berry_arr>0);
    spikePos = spikePos(berr_log);
    
    for spikeCount = 1:length(spikePos)
        ber = berry_arr(spikePos(spikeCount));
        
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4],c,'LineWidth',...
            ber*0.1);
    end
end

ylim([0 size(berry_cell, 1)+1]);