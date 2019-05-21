function [seg_cortex,seg_sub] = F_dual_segment(fun_data,mask_cortex,mask_sub,seg_initial,seg_thre,seg_strategy)
%dual_segment
% dual_segment of subcortical regions
% fun_data: 4-D functional fMRI matrix
% mask_cortex: mask of cortex regions
% mask_sub: mask of subcortical regions
% seg_initial: input the number of sub regions or initial segment of cortex
% seg_thre: the threshold to terminate the segment
% seg_strategy: 1:group 2:individual
% region or subcortical region
%
% created by Heng, 2018/10/16

if(nargin == 5)
    seg_strategy = 1;
end

if(iscell(fun_data) && seg_strategy == 2)
    seg_cortex = {};
    seg_sub = {};
    
    for i = 1:length(fun_data)
        disp(['processing: ',num2str(i),'/',num2str(length(fun_data))]);
        [seg_cortex_tmp,seg_sub_tmp] = F_dual_segment(fun_data{i},...
            mask_cortex,mask_sub,seg_initial,seg_thre);
        seg_cortex = cat(1,seg_cortex,seg_cortex_tmp);
        seg_sub = cat(1,seg_sub,seg_sub_tmp);
    end
    return;
end

if(isscalar(seg_initial))
    seg_num = seg_initial;
    seg_initial_cortex = mask_cortex;
    seg_initial_cortex(seg_initial_cortex>0) = ...
        randi(seg_num,1,sum(seg_initial_cortex(:)>0));
elseif(all((seg_initial(:)>0) == (mask_cortex(:)>0)))
    seg_num = max(seg_initial(:));
    seg_initial_cortex = seg_initial;
elseif(all((seg_initial(:)>0) == (mask_sub(:)>0)))
    seg_num = max(seg_initial(:));
    seg_initial_sub = seg_initial;
else
    error('initial segment does not match masks!!');
end

if(~exist('seg_initial_cortex','var'))
    seg_initial_cortex = F_cluster(fun_data,seg_initial_sub,mask_cortex,seg_num);
end

[seg_cortex,seg_sub] = F_dual_cluster(fun_data,seg_initial_cortex,mask_sub,seg_num,seg_thre);
disp('done');
end

function [seg_cortex_i,seg_sub_i] = F_dual_cluster(fun_data,seg_cortex_i_b,mask_sub,seg_num,seg_thre)
mask_cortex_i = double(seg_cortex_i_b>0);
num_cortex = sum(mask_cortex_i(:)>0);
for i_1 = 1:10
    seg_sub_i = F_cluster(fun_data,seg_cortex_i_b,mask_sub,seg_num);
    seg_cortex_i_a = F_cluster(fun_data,seg_sub_i,mask_cortex_i,seg_num);
    match_rate = ...
        sum((seg_cortex_i_a(:) == seg_cortex_i_b(:)) & (mask_cortex_i(:)>0))./num_cortex;
    if(match_rate >=seg_thre)
        seg_cortex_i = seg_cortex_i_a;
        return;
    else
        disp(['loop: ',num2str(i_1),...
            ' the similarity with previous loop is ',num2str(match_rate)]);
        seg_cortex_i_b = seg_cortex_i_a;
        seg_cortex_i = seg_cortex_i_b;
    end
end
warning('match rate did not reach 100% after 100 loops');
end

function seg_b = F_cluster(fun_data,seg_a,mask_b,seg_num)
    
    if(iscell(fun_data))
        for i_2 = 1:length(fun_data)
            fun_data_1d_2 = reshape(fun_data{i_2},...
                size(fun_data{i_2},1)*size(fun_data{i_2},2)*size(fun_data{i_2},3),...
                size(fun_data{i_2},4));
            seg_a_1d = reshape(seg_a,length(seg_a(:)),1);
            seg_num_2 = max(seg_a_1d);
            for i_seg_num_2 = 1:seg_num_2
                a_sig(i_seg_num_2,:) = mean(fun_data_1d_2(seg_a_1d == i_seg_num_2,:),1);
            end
            mask_b_1d = reshape(mask_b,length(mask_b(:)),1);
            if(i_2 == 1)
                connect_matrix_2 = corr(a_sig',fun_data_1d_2(mask_b_1d>0,:)');
            else
                connect_matrix_2 = connect_matrix_2 + corr(a_sig',fun_data_1d_2(mask_b_1d>0,:)');
            end
            clear a_sig
        end
    else
        fun_data_1d_2 = reshape(fun_data,...
            size(fun_data,1)*size(fun_data,2)*size(fun_data,3),size(fun_data,4));
        seg_a_1d = reshape(seg_a,length(seg_a(:)),1);
        seg_num_2 = max(seg_a_1d);
        for i_2 = 1:seg_num_2
            a_sig(i_2,:) = mean(fun_data_1d_2(seg_a_1d == i_2,:),1);
        end
        mask_b_1d = reshape(mask_b,length(mask_b(:)),1);
        connect_matrix_2 = corr(a_sig',fun_data_1d_2(mask_b_1d>0,:)');
        clear a_sig
    end
    [~,seg_b_tmp] = max(connect_matrix_2);
    for i_2 = 1:seg_num
        voxel_num_i(i_2) = sum(seg_b_tmp(:) == i_2);
    end
%     if((max(voxel_num_i)./min(voxel_num_i)) > 10)
%         [~,idx_max] = max(voxel_num_i);
%         [~,idx_min] = min(voxel_num_i);
%         seg_b_tmp_max = seg_b_tmp(seg_b_tmp == idx_max);
%         seg_b_tmp_max(rand(size(seg_b_tmp_max))>0.5) = idx_min;
%         seg_b_tmp(seg_b_tmp == idx_max) = seg_b_tmp_max;
%     end
    seg_b_1d = zeros(size(mask_b_1d));
    seg_b_1d(mask_b_1d>0) = seg_b_tmp;
    seg_b = reshape(seg_b_1d,size(mask_b)); 
end

    