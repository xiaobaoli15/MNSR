%Embedding_space
function [Gallery,Probe]=Embedding_space(Training_data,Training_labels,means_Gallery,means_Probe,embedding_space_method,a_train_Labels,b_train_Labels,means_a_Train,means_b_Train,cluster_num,numClasses)
     switch embedding_space_method
        case 'Unsupervised'
             fprintf('%s内嵌空间----->>',embedding_space_method);
             [Gallery,~]=normalizeBase(means_Gallery);
             [Probe,~]=normalizeBase(means_Probe);
        case 'KLFDA'
             fprintf('%s内嵌空间----->>',embedding_space_method);
             AlgoOption.kernel ='linear';  %chi2,exp,poly,chi2-rbf,linear                                                                                                                                                                                                                                                                                                                     
             AlgoOption.npratio =0; % npratio is not required.
             AlgoOption.beta =0.0001;%0.0001
             AlgoOption.d =500;
             AlgoOption.epsilon =1e-6;
             AlgoOption.LocalScalingNeighbor =8;%8
             [Method,~]= kLFDA(Training_data',Training_labels',AlgoOption);%LFDA算法效果更好一点
             [ker_Gallery] = ComputeKernelTest(Training_data',means_Gallery', Method);
             [ker_Probe] = ComputeKernelTest(Training_data',means_Probe', Method);
             T=Method.P;
             [Gallery,~]=normalizeBase(T*ker_Gallery);
             [Probe,~]=normalizeBase(T*ker_Probe);
        case 'KISSME'
             fprintf('%s内嵌空间----->>',embedding_space_method);
             options.N=size(a_train_Labels,2); % not used actually
             options.lambda=0.001;
             options.npratio = 3;%正负样本的比例
             cHandle=LearnAlgoKISSME(options);
             [pos_pairs,negtive_pairs]=split_neg_pos(numClasses,cluster_num);
             pairs=[pos_pairs;negtive_pairs];
             y = [ones(size(a_train_Labels,2),1); ones(size(b_train_Labels,2),1).*(-1)];
             s=learnPairwise(cHandle,double(Training_data),pairs(:,1),pairs(:,2),y>0);
             ds.(cHandle.type)=s;
             T=ds.kissme.M;
             L=chol(T);
             [Gallery,~]=normalizeBase(L*means_Gallery);
             [Probe,~]=normalizeBase(L*means_Probe);                   
        case 'XQDA'
             fprintf('%s内嵌空间----->>',embedding_space_method);
             [W, M] = XQDA(means_a_Train' ,means_b_Train',a_train_Labels',b_train_Labels');
             L=chol(M);
             [Gallery,~]=normalizeBase(L*W'*means_Gallery);
             [Probe,~]=normalizeBase(L*W'*means_Probe);
        case 'MFA'
             fprintf('%s内嵌空间----->>',embedding_space_method);
             AlgoOption.Nw = 0; % 0--use all within class samples
             AlgoOption.Nb = 12;
             AlgoOption.d = 100;
             AlgoOption.beta = 0.01;
             AlgoOption.epsilon =1e-6;
             AlgoOption.kernel ='linear';  %chi2,linear
             AlgoOption.dataname ='';  %chi2  
             [Method,~]=  MFA(Training_data',Training_labels',AlgoOption);
             [ker_Gallery] = ComputeKernelTest(Training_data',means_Gallery', Method);
             [ker_Probe] = ComputeKernelTest(Training_data',means_Probe', Method);
             T=Method.P;
             [Gallery,~]=normalizeBase(T*ker_Gallery);
             [Probe,~]=normalizeBase(T*ker_Probe);     
        case 'LMNN'
            options.npratio = 3;
            cHandle=LearnAlgoLMNN();
            [pos_pairs,negtive_pairs]=split_neg_pos(numClasses,cluster_num);
            pairs=[pos_pairs;negtive_pairs];
            y = [ones(size(a_train_Labels,2),1); ones(size(b_train_Labels,2),1).*(-1)];
            s=learnPairwise(cHandle,Training_data,pairs(:,1),pairs(:,2),y>0);
            T=ds.lmnn.M;  
            [Gallery,~]=normalizeBase(T*means_Gallery);
            [Probe,~]=normalizeBase(T*means_Probe);     
         case 'ITML'
            options.npratio = 3;
            cHandle=LearnAlgoITML();
            [pos_pairs,negtive_pairs]=split_neg_pos(numClasses,cluster_num);
            pairs=[pos_pairs;negtive_pairs];
            y = [ones(size(a_train_Labels,2),1); ones(size(b_train_Labels,2),1).*(-1)];
            s=learnPairwise(cHandle,Training_data,pairs(:,1),pairs(:,2),y>0);
            ds.(cHandle.type)=s;
            T=ds.itml.M;
            [Gallery,~]=normalizeBase(T*means_Gallery);
            [Probe,~]=normalizeBase(T*means_Probe);    
        otherwise
            fprintf('没有该方法，请从Usupervised,KLFDA,KISSME,XQDA,MFA,ITML,LMNN 中选择一种。');
     end
end