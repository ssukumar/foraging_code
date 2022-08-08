function run_all_subjs ( exp_no)

if nargin<1
    exp_no=1;
end

filename = 'subj_order.txt';
folder = '..';
filename = fullfile(folder, filename);
subj_order = readtable(filename, 'Delimiter',' ',...
    'ReadVariableNames',0);
subj_order = table2cell(subj_order);

if exp_no ==1
    subjs = {subj_order{:,2}};
    n_subjs = length(subjs);
    for s  = 1:n_subjs 
        
        subj = subjs{s};
        analyze_foraging_data(subj,1,1);
        
    end

elseif exp_no == 2
    subjs = {'XGZ','EDS','LWO',...
    'SHD','BYT','ALT','CKL', 'PJS','NTK'};
    n_subjs = length(subjs);
    for s  = 1:n_subjs 
        
        subj = subjs{s};
        analyze_foraging_data(subj,exp_no,1);
        
    end
    
else 
     subjs = {'ZME','NSY','XGZ','AVX','EDS','LWO',...
    'SHD','BYT','ALT','CKL', 'PJS','NTK'};
    n_subjs = length(subjs);
    for s  = 1:n_subjs 
        
        subj = subjs{s};
        analyze_fml_foraging(subj);
        
    end
end

