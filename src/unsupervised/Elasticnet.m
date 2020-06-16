%Elasticnet
function [Alphas,E]=Elasticnet(Gallery,Probe)
    lambda1 = 0.1;
    lambda2 =0.1;
    opts.loss = 'l2'; 
    [Alphas,E,~,~,~] = elasticnetR(Gallery,Probe,lambda1,lambda2,opts);
end