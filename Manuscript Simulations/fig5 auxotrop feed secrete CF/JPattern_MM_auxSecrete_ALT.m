function M = JPattern_MM_auxSecrete_ALT(ncells,n,p,cell_type)
	% Number of auxotrophies in one cell = 1
	p1 = p-1;
	each_cell = (n*p1+1)+1; 	
	tot_siz = each_cell*ncells;  % Size without crossfeed terms
    first_S = 1:each_cell:tot_siz;
	
	% With Secretion
	M = zeros(tot_siz + 1);	% +1 eq for the total final secreted metabolite
		
	for ncell = 1:ncells
		% Main Enzyme reactions
		index = first_S(ncell)+(0:(n-1)*p1);
		M(index,index) = eye(length(index));
		
		% Common Substrate Dependence (of S2A, S2B...on S1)
		index = first_S(ncell)+(1:p1);
		M(index,first_S(ncell)) = ones(p1,1);
		
		% Intermediate cascade reactions (of S3A, S3B... depending upon S2A, S2B...)
		index = first_S(ncell)+p1+(1:(n-2)*p1);
		M(index, index-p1) = M(index, index-p1) + diag(ones((n-2)*p1,1));
	
		sec_n = cell_type(2);
		sec_p = cell_type(3);
		aux_p = cell_type(4);
		
		% Final Metabolites may only be produced or uptaken. Secretion of final metabolite is instead considered from the previous step metabolite
		index = first_S(ncell) + (n-1)*p1 + setdiff(1:p,aux_p);
		M(index, index-p1) = M(index, index-p1) + diag(ones(p1,1));
		
		% Secreted metabolite external amount variable
		M(tot_siz + 1, first_S(ncell) + ((sec_n-1)-1)*p1 + (sec_p>p1)*(p1-sec_p)+sec_p) = 1;
		
	end
	
%	M = sparse(M);
%	display(M);
end
		