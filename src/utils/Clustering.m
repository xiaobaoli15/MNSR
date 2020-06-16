%% clustering
%each_sequence_person_num 
%Data                         :requried data
%cluster_num                  :the number of clusters we need
%numClasses                   :the number of identities
%do_cluster                   :a flag indacates whether to cluster or not
%%
function clustered=Clustering(Data,cluster_num,numClasses,each_sequence_person_num,do_cluster)
     clustered=[];
     cur_sum=1; 
     switch do_cluster
         case '0'
             %fprintf('Without clustering----->>')
             for j=1:numClasses
                 shuffle_num = randperm(each_sequence_person_num(j));
                 img_selected = shuffle_num(1:cluster_num);
                 cur_sequence=Data(:, cur_sum+img_selected-1);
                 clustered=cat(2,clustered,cur_sequence);
                 cur_sum= cur_sum+ each_sequence_person_num(j);
             end
         case '1'
             %fprintf('K-means clustering----->>')
             for j=1:numClasses  
                 cur_sequence=Data(:,cur_sum:cur_sum+each_sequence_person_num(j)-1);
                 [~,MEANS] = kmeans(cur_sequence',cluster_num);
                 clustered=cat(2,clustered,MEANS');
                 cur_sum=cur_sum+each_sequence_person_num(j);   
             end
         otherwise
            fprintf('the clustering parameter error:you can set do_cluster as 0 or 1.');
     end
end
