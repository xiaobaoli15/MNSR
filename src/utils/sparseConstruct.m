function sparse_rapresentation_res=sparseConstruct(Probe,Gallery,Alphas,idPersons,Gallery_Labels,E)
    errors=[];
    for t=1:size(Probe,2)
        Alpha=Alphas(:,t);
        each_Probe=Probe(:,t);
        eachE=E(:,t);
        each_Probe=each_Probe-eachE;
        Alphaband=repmat(Alpha,1,length(idPersons));
        for p=1:length(idPersons)
            idx=find(idPersons(p)~=Gallery_Labels & Gallery_Labels~=0);
            Alphaband(idx,p)=0;
        end
        residual=computeResidualColumnWise(each_Probe,idPersons,Gallery,Alphaband); 
        error=normColumnWise(residual,2);
        error=normalizeBase(error);  
        errors=cat(2,errors,error);
    end
    sparse_rapresentation_res=errors;
end
