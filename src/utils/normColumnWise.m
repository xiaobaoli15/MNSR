%normColumnWiseº¯Êý
function n=normColumnWise(residual,p)
       n=sum(residual.^p,1).^(1/p);
       n=n';
return