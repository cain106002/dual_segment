function seg_num = F_determine_seg_num(fun_data,mask_sub,mask_cortex)


mask_sub_1d = reshape(mask_sub,...
    length(mask_sub(:)),1);
mask_cortex_1d = reshape(mask_cortex,...
    length(mask_cortex(:)),1);
if(iscell(fun_data))
    corr_matrix = zeros(sum(mask_sub_1d>0),sum(mask_cortex_1d>0));
    for i = 1:length(fun_data)
        fun_data_1d = reshape(fun_data{i},...
        size(fun_data{i},1)*size(fun_data{i},2)*size(fun_data{i},3),size(fun_data{i},4));
        sub_sig = fun_data_1d(mask_sub_1d>0,:);
        cortex_sig = fun_data_1d(mask_cortex_1d>0,:);
        corr_matrix = corr_matrix+corr(sub_sig',cortex_sig');
    end
    corr_matrix = corr_matrix./length(fun_data);
else
    fun_data_1d = reshape(fun_data,...
        size(fun_data,1)*size(fun_data,2)*size(fun_data,3),size(fun_data,4));
    sub_sig = fun_data_1d(mask_sub_1d>0,:);
    cortex_sig = fun_data_1d(mask_cortex_1d>0,:);
    corr_matrix = corr(sub_sig',cortex_sig');
end
eva = evalclusters(corr_matrix,'kmeans','Silhouette','kList',[2:20],...
    'Distance','correlation');
seg_num = eva.OptimalK;
figure
plot(eva.InspectedK,eva.CriterionValues);
end