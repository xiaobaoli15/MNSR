%for normalize
function [Dnorm,norms]=normalizeBase(D)
    norms=sqrt(sum(D.^2));
    norm_band=repmat(norms,size(D,1),1);
    Dnorm=D./norm_band;
    %Dnorm=mexNormalize(D);
return
