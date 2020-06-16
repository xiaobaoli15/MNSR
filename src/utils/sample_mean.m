function dist=sample_mean(sample_mat1,sample_mat2,cluster_num) 
    mean_mat1=[];
    mean_mat2=[];
    %mean
    for i=1:size(sample_mat1,2)/cluster_num
        each_avg_mat1=mean(sample_mat1(:,(i-1)*cluster_num+1:i*cluster_num),2);
        each_avg_mat2=mean(sample_mat2(:,(i-1)*cluster_num+1:i*cluster_num),2);
        mean_mat1=[mean_mat1 each_avg_mat1];
        mean_mat2=[mean_mat2 each_avg_mat2];
    end
    %dist
    dist=[];
    for i=1:size(mean_mat1,2)
        diff=mean_mat1(:,i)-mean_mat2(:,i);
        tmp_dist=sqrt(sum(diff.^2));
        dist=[dist tmp_dist];
    end