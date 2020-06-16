  %2018.7.11
  %xiaobaili15@163.com
  %bi-sparse
  function [avg_Alphas]=bi_sparse(Gallery,Probe,cluster_num)	  
	  [Alphas,E]=L1_norm(Gallery,Probe);
	  avg_Alphas=[];
	  for i=1:size(Probe,2)/cluster_num
		  Alphas_band=Alphas(:,(i-1)*cluster_num+1:i*cluster_num);
		  avg_Alphas_band=means(Alphas_band,2);
		  avg_Alphas=[avg_Alphas avg_Alphas_band];
	  end
  end