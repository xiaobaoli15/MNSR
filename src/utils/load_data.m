%load data
if strcmp(Dataset,'PRID 2011')
    fprintf('Dataset:%s --->>',Dataset);
    load('features/PRID 2011/a_each_num.mat');
    load('features/PRID 2011/b_each_num.mat');
    load('database/train-test people splits/train_test_splits_prid.mat');%the split protocal for repeate 10 trails
    switch feature_method
        case 'HistLBP'
             fprintf('Feature:%s----->>',feature_method);
             load('features/PRID 2011/HistLBP/test_features.mat');
             load('features/PRID 2011/HistLBP/train_features.mat');
             a_features=test_feature';%2580*43477
             b_features=train_feature';%2580*51510
        case 'LOMO'
             fprintf('Feature:%s----->>',feature_method);
             load('features/PRID 2011/LOMO/test_descriptors_LOMO.mat');
             load('features/PRID 2011/LOMO/train_descriptors_LOMO.mat');
             a_features=test_descriptors;%26960*43477
             b_features=train_descriptors;%26960*51510
        case 'GOG'
             fprintf('Feature:%s----->>',feature_method);
             load('features/PRID 2011/GOG/test_GOG_features_7.mat');
             load('features/PRID 2011/GOG/train_GOG_features_7.mat');
             a_features=test_GOG';%16215*43477
             b_features=train_GOG';%16215*51510
        otherwise
            fprintf('Without this discriptor,please choose one from HistLBP,LOMO,GOG');
    end
elseif (strcmp(Dataset,'iLIDS-VID'))
    fprintf('Dataset:%s --->>',Dataset);
    load('features/iLIDS-VID/a_each_num.mat');%1*300,a_each_num
    load('features/iLIDS-VID/b_each_num.mat');%1*300,b_each_num
    load('database/train-test people splits/train_test_splits_ilidsvid.mat');%ls_set,10*300
    switch feature_method
        case 'HistLBP'
             fprintf('Feature:%s----->>',feature_method);
             load('features/iLIDS-VID/HistLBP/test_features.mat');
             load('features/iLIDS-VID/HistLBP/train_features.mat');
             a_features=test_feature';%2580*19701
             b_features=train_feature';%2580*22758
        case 'LOMO'
             fprintf('Feature:%s----->>',feature_method);
             load('features/iLIDS-VID/LOMO/test_descriptors_LOMO.mat');%26960*19701
             load('features/iLIDS-VID/LOMO/train_descriptors_LOMO.mat');%26960*22758
             a_features=test_descriptors;%26960*19701
             b_features=train_descriptors;%26960*22758
        case 'GOG'
             fprintf('Feature:%s----->>',feature_method);
             load('features/iLIDS-VID/GOG/test_GOG_features_7.mat');%19701*7567
             load('features/iLIDS-VID/GOG/train_GOG_features_7.mat');%22758*7567
             a_features=test_GOG';%7567*19701
             b_features=train_GOG';%7567*22758
        otherwise
           fprintf('Without this discriptor,please choose one from HistLBP,LOMO,GOG');
    end
elseif (strcmp(Dataset,'SAIVT-SoftBio'))
    fprintf('Dataset:%s --->>',Dataset);
    load('features/SAIVT-SoftBio/camA_each_sequence_num.mat');
    load('features/SAIVT-SoftBio/camB_each_sequence_num.mat');
    a_each_num=camA_each_sequence_num;
    b_each_num=camB_each_sequence_num;
    load('database/train-test people splits/train_test_splits_saivt_softbio.mat');
    switch feature_method
        case 'HistLBP'
             fprintf('Feature:%s----->>',feature_method);
             load('features/SAIVT-SoftBio/HistLBP/test_features.mat');
             load('features/SAIVT-SoftBio/HistLBP/train_features.mat');
             a_features=test_feature';%2580*43477
             b_features=train_feature';%2580*51510
        case 'LOMO'
             fprintf('Feature:%s----->>',feature_method);
             load('features/SAIVT-SoftBio/LOMO/test_descriptors_LOMO.mat');
             load('features/SAIVT-SoftBio/LOMO/train_descriptors_LOMO.mat');
             a_features=test_descriptors;
             b_features=train_descriptors;
        case 'GOG'
             fprintf('Feature:%s----->>',feature_method);
             load('features/SAIVT-SoftBio/GOG/camA_feature.mat');
             load('features/SAIVT-SoftBio/GOG/camB_feature.mat');
             a_features=camA_GOG_features;
             b_features=camB_GOG_features;
        otherwise
           fprintf('Without this discriptor,please choose one from HistLBP,LOMO,GOG');
    end
else
    fprintf('Without this dataset,please choose one from PRID 2011,iLIDS-VID,SAIVT-SoftBio.')
end
        