%KISSME算法的正负样本对的划分为i:
function [pos_pairs,negtive_pairs]=split_neg_pos(g_num,cluster_num)
N=g_num*cluster_num;%数据集中总的行人图像数量
N_sequence=1:N;%图像标签
negtive_pairs=[];
%产生负样本对
for i=1:g_num
    relative_positive=[(i-1)*cluster_num+1:i*cluster_num];
    relative_negtive=setdiff(N_sequence,relative_positive);%返回在A中有，而B中没有的值，结果向量将以升序排序返回
    neg_sequence=N+relative_negtive;
    rand_sequence_neg=neg_sequence(randperm(length(neg_sequence)));
    neg=rand_sequence_neg(1:cluster_num);
    rand_sequence_pos=relative_positive(randperm(length(relative_positive)));
    pos=rand_sequence_pos(1:cluster_num);
    tmp_negtive_pair=[pos' neg'];
    negtive_pairs=[negtive_pairs;tmp_negtive_pair];
end

%产生正样本对
a_index=1:N;
b_index=N+1:2*N;
tmp_pos_pairs=[a_index;b_index];
pos_pairs=tmp_pos_pairs';
end
    
    

    