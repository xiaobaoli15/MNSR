%for compute acc of label estimation
function label_estimation_acc(a1,a2,estimate_a,estimate_b)
     count=0;
     for i=1:length(estimate_a)
         index1=find(a1(estimate_a(i))==a2);
         if index1==estimate_b(i)            
             count=count+1;
         end
     end
     fprintf('')
     fprintf('Total num:[%d]|acc num:[%d]--->>',length(estimate_a),count)    
end