function M = JPattern_noisy_cat(ncells,n,p,cf_enable)
	
	tot_siz = (p*n+1)*ncells;  % Size without crossfeed terms
    first_S = 1:(n*p+1):tot_siz;
	
	if cf_enable == 0
		% No CrossFeeding
		M = zeros(tot_siz,tot_siz);
	elseif cf_enable == 1
		% With CrossFeeding
		M = zeros(tot_siz+n*p,tot_siz+n*p);
		
		% CrossFeeding Self Dependence - Uptake
		index = tot_siz+(1:n*p);
		M(index, index) = eye(n*p);
	elseif cf_enable == 2
		% With CrossFeeding and Analysis enabled
		M = zeros(tot_siz+3*n*p, tot_siz+3*n*p);
		
		% CrossFeeding Self Dependence - Uptake
		index = tot_siz+(1:n*p);
		M(index, index) = eye(n*p);
		
		% ANALYSIS CrossFeeding Self Dependence - Uptake
		M(index+2*n*p, index) = eye(n*p);
	end
	
	% Main Enzyme reactions
	index = 1:tot_siz;
	M(index,index) = eye(tot_siz);

	for cell = 1:ncells
		% Multiple cascade reactions
		index = first_S(cell)+(0:n*p);
		M(index, index) = M(index, index) + diag(ones(((n-1)*p+1),1),-p);
		
		% Common Substrate Dependence 
		index = first_S(cell)+(1:p);
		M(index,first_S(cell)) = ones(p,1);
		
		if cf_enable > 0
			% CrossFeeding Steps (Rows)
			index = first_S(cell)+(1:p*n);
			M(index, tot_siz+(1:n*p)) = M(index, tot_siz+(1:n*p)) + eye(n*p);
			
			% CrossFeeding Steps (Columbs)
			index = first_S(cell)+[1:p*n];
			M(tot_siz+(1:n*p), index) = M(tot_siz+(1:n*p), index) + eye(n*p);
			
			if cf_enable == 2
				% ANALYSIS: CrossFeeding Steps (Columbs)
				M(tot_siz+n*p+(1:n*p), index) = M(tot_siz+n*p+(1:n*p), index) + eye(n*p);
			end
		end
	end
	% Add the extra variable of noisy catabolite feed
	M = blkdiag(1,M);
	M(2,1) = 1;
	M = sparse(M);
%	display(M);
end