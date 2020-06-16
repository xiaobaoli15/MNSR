function [a_label,b_label]=label_estimation(X1,X2,cluster_num,g_num,K,num_selected)
    probablity1=generate_avg_probablity(X1,cluster_num,g_num);
    probablity2=generate_avg_probablity(X2,cluster_num,g_num);
    
    %sort
    [value1,idx1]=sort(probablity1,'descend');
    [value2,idx2]=sort(probablity2,'descend');
    
    %KNN
    KNN_idx1=idx1(1:K,:);
    KNN_value1=value1(1:K,:);
    KNN_idx2=idx2(1:K,:);
    KNN_value2=value2(1:K,:);
    
    %mutual KNN
    a_label=[];
    b_label=[];
    mutual_value_list=[];
    [m,n]=size(probablity1);
    %a->b
    for i=1:n
        mutual_KNN_idx_list=[];
        mutual_KNN_value_list=[];
        each_KNN_idx1=KNN_idx1(:,i);
        each_KNN_value1=KNN_value1(:,i);
        for j=1:K
            each_KNN_idx2=KNN_idx2(:,each_KNN_idx1(j));
            each_KNN_value2=KNN_value2(:,each_KNN_idx1(j));
            if ismember(i,each_KNN_idx2)
               mutual_KNN_idx_list=[mutual_KNN_idx_list each_KNN_idx1(j)]; 
               mutual_KNN_value_list=[mutual_KNN_value_list (each_KNN_value2(find(i==each_KNN_idx2))+each_KNN_value1(j))/2];
            end
        end

        if  size(mutual_KNN_idx_list)==0
            continue;
        else
            [mutual_value,mutual_idx]=max(mutual_KNN_value_list);
            selected_idx=mutual_KNN_idx_list(mutual_idx); 
            
            mutual_value_list=[mutual_value_list mutual_value];
            
            a_label=[a_label i];
            b_label=[b_label selected_idx];
        end      
    end
    [~,idx1]=sort(mutual_value_list,'descend');
    if length(mutual_value_list)<=num_selected
       num_selected=length(mutual_value_list);        
    end
        
    a_label=a_label(idx1(1:num_selected));
    b_label=b_label(idx1(1:num_selected));
    
end
