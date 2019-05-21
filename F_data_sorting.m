function data_sorted = F_data_sorting(data,data_ref,seg_num)
idx_perm = perms([1:seg_num]);
voxel_num = sum(data(:)>0);
match_rate_max = 0;
for i = 1:length(idx_perm)
    data_tmp = zeros(size(data));
    for j = 1:seg_num
        data_tmp(data == j) = idx_perm(i,j);
    end
    match_rate = zeros(length(data_ref),1);
    for j = 1:length(data_ref)
        match_rate(j) = sum(data_tmp(data_tmp>0) == data_ref{j}(data_tmp>0))./voxel_num;
    end
    match_rate_mean = mean(match_rate);
    if(match_rate_mean>match_rate_max)
        match_rate_max = match_rate_mean;
        data_sorted = data_tmp;
    end
end