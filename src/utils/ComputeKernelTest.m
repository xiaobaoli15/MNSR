function [K_test] = ComputeKernelTest(train, test, Method)
if (size(train,2))>2e4 && (strcmp(Method.kernel, 'chi2') || strcmp(Method.kernel, 'chi2-rbf'))
    % if the input data matrix is too large then use parallel computing
    % tool box.
%     poolobj = parpool;
    
    switch Method.kernel
        case {'linear'}
            K_test = train * test';
        case {'chi2'}
            for i =1:size(test,1)
                dotp = bsxfun(@times, test(i,:), train);
                sump = bsxfun(@plus, test(i,:), train);
                K_test(:,i) = 2* sum(dotp./(sump+1e-10),2);
            end
        case {'chi2-rbf'}
            sigma = Method.rbf_sigma;
            for i =1:size(test,1)
                subp = bsxfun(@minus, test(i,:), train);
                subp = subp.^2;
                sump = bsxfun(@plus, test(i,:), train);
                K_test(:,i) =  sum(subp./(sump+1e-10),2);
                %disp(i);
            end
            K_test =exp(-K_test./sigma);
    end
%     delete(poolobj)
else
    switch Method.kernel
        case {'linear'}
            K_test = train * test';
        case {'chi2'}
            for i =1:size(test,1)
                dotp = bsxfun(@times, test(i,:), train);
                sump = bsxfun(@plus, test(i,:), train);
                K_test(:,i) = 2* sum(dotp./(sump+1e-10),2);
               % disp(i);
            end
        case {'chi2-rbf'}
            sigma = Method.rbf_sigma;
            for i =1:size(test,1)
                subp = bsxfun(@minus, test(i,:), train);
                subp = subp.^2;
                sump = bsxfun(@plus, test(i,:), train);
                K_test(:,i) =  sum(subp./(sump+1e-10),2);
               % disp(i);
            end
            K_test =exp(-K_test./sigma);
    end
end
return;