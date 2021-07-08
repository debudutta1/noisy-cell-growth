% Function to create a Matrix form ODE for simulation of growth of auxotrophic cells crossfeeding with each other
function dx = ODE_gen_MM_auxFeedSecrete_ALT(t,x, ncells, p, n, Sin, par, metab_threshold, cell_dat, n_ext)
	
	% MATLAB supposedly passes by reference by default if the variable is never modified in the function!
	ET = zeros(ncells,p,n);
	% Interpolate the value of Total Enzyme at given time
	for cell = 1:ncells
		for parcas = 1:p
			for sercas = 1:n
				ET(cell,parcas,sercas) = interpl_lastval(cell_dat{cell}.prot_dat{parcas,sercas}.cum_prot_profile.t,cell_dat{cell}.prot_dat{parcas,sercas}.cum_prot_profile.v,t);
			end
		end
	end
	
	neqs = (p*n+1)+1;	 	% +1 Since uptake missing metab
	nterms = (p*n+1);		% No need for +1 Since secretion is now a direct ratio of production
	
	% Set up total Matrix and Vector
	M = zeros(neqs*ncells + 1, nterms*ncells);	% +1 eq for the total final secreted metabolite
	V = zeros(nterms*ncells, 1);	% Secretion doesn't need more global terms, already in nterms. Uptake is constant
	
	first_S = 1:neqs:ncells*neqs;	% Equation Index of the first Substrate of each cell
	first_V = 1:nterms:ncells*nterms;	% First Term Index of each cell in V
	
	%% M
	count = 0;
	for ncell = 1:ncells
		% Cell type
		c_type = cell_dat{ncell}.type(1);
		
		%% M
		M(first_S(ncell),first_V(ncell)) = Sin;
		
		% Auxotrophy direct import
		M(first_S(ncell)+n*p+1, first_V(ncell)) = par.feed_rate; % Since Missing metabolite is placed as the last term always! first_V(ncell) since directly imported without any additional dependent variable
		
		% Production terms in M. all cells are single auxotrophs (p = p0-1)
		len1 = first_S(ncell)+(1:n*p);
		len2 = first_V(ncell)+(1:n*p);
		M(len1,len2) = M(len1,len2) + diag(ET(ncell,:).*par.k3(c_type,:));	
		
		% Consumption terms in M
		M(first_S(ncell),first_V(ncell)+(1:p)) = ET(ncell,:,1).*-par.k3(c_type,:,1);
		if n > 1
			len1 = [first_S(ncell)+(1:p*(n-1))];
			len2 = [first_V(ncell)+p+(1:p*(n-1))];
			M(len1,len2) = M(len1,len2) + diag(reshape(ET(ncell,:,2:n),1,[]).*reshape(-par.k3(c_type,:,2:n),1,[]));
		end
		
		% Decide what to secrete based on cell type
		% [1, 2, 3, 4]. 1= cell type, (2,3) = (n,p) index of secreted metab, 4 = p of auxoptrophy 
		% Index of secreted and auxotrophy
		sec_n = cell_dat{ncell}.type(2);
		sec_p = cell_dat{ncell}.type(3);
		aux_p = cell_dat{ncell}.type(4);
		
		% Per Cell Secretions
		% Adjust final metab and Per Cell secretion
		% Case 1 - Secrete Last Metab in cascade % sec_n = n  % Case 2 - Secrete Second Last Metab in cascade, sec_n = n-1
		len1 = first_S(ncell) + (sec_n-1)*p + (sec_p>p)*(p-sec_p)+sec_p;
		temp_var = M(len1,len1);
		M(len1,len1) = M(len1,len1)*(1-par.secRatio);	% Reduce metab by (1-secRatio) fraction. Add secRatio to secretion
		
		% Cell Secretions Added to External [Metabolite] 
		len1 = ncells*neqs + c_type;				% Ext Metab Eq
		len2 = first_S(ncell) + (sec_n-1)*p + (sec_p>p)*(p-sec_p)+sec_p;
		M(len1,len2) = M(len1,len2) + temp_var*par.secRatio;		
		
		
		%% V
		V(first_V(ncell)) = 1;
		
		% Common Substrate step
		V(first_V(ncell)+(1:p)) = x(first_S(ncell))./( x(first_S(ncell)) + par.KM(1:p,1) );
		
		% Subsequent Steps in cascade
		if n > 1
			V(first_V(ncell)+p+(1:p*(n-1))) = x(first_S(ncell)+(1:p*(n-1)))./( x(first_S(ncell)+(1:p*(n-1))) + reshape(par.KM(:,2:n),[],1) );
		end
		
	end		
	
% 	V' 
% 	M
	dx = M*V;
end