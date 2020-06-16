function [postive_value,postive_idx]=mean_set_dist(Gallery,Probe,g_num,cluster_num)
    probe_gallery_dist=[];
    for euc_dist=1:size(Probe,2)
        minu=bsxfun(@minus,Probe(:,euc_dist),Gallery);
        probe_gallery_dist=[probe_gallery_dist;sqrt(sum(minu.^2,1))];       
    end

    %average the distance of  probe images with same identites
    col_avg_dist=[];
    row_avg_dist=[];
    for j=1:g_num
            each_avg_X=sum(probe_gallery_dist(:,(j-1)*cluster_num+1:j*cluster_num),2);
            col_avg_dist=[col_avg_dist each_avg_X];
    end
    for i=1:g_num
        each_avg_X=sum(col_avg_dist((i-1)*cluster_num+1:i*cluster_num,:),1);
        row_avg_dist=[row_avg_dist;each_avg_X];    
    end
    mean_dist=row_avg_dist;
    mean_dist=mean_dist';%每一列是一个probe行人图像序列
    
    %sort mean_dist
    [~,idx1]=sort(mean_dist,'ascend');
    postive_value=1:size(Probe,2);
    postive_idx=idx1(1,:);
end