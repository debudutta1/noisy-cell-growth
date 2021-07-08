% Function to create a Matrix form ODE for simulation of growth of cells in terms of metabolite production only. NO CROSS-FEEDING
function dx = ODE_gen_MMfn_analysis_extFeed(t,x, ncells, p, n, Sin, par, metab_threshold, cell_dat, n_ext)
	% MATLAB supposedly passes by reference by default if the variable is never modified in the function!
	ET = zeros(ncells,p,n);
	% Interpolate the value of Total Enzyme at given time
	for ncell = 1:ncells
		for parcas = 1:p
			for sercas = 1:n
				%ET(ncell,parcas,sercas) = interpl_lastval(cell_dat{ncell}.prot_dat{parcas,sercas}.tot_prot_profile.t,cell_dat{ncell}.prot_dat{parcas,sercas}.tot_prot_profile.v,t);
				ET(ncell,parcas,sercas) = interpl_lastval(cell_dat{ncell}.prot_dat{parcas,sercas}.cum_prot_profile.t,cell_dat{ncell}.prot_dat{parcas,sercas}.cum_prot_profile.v,t);
			end
		end
	end
% 	ncells = 2; n = 3; p = 3;
% 	ET = 1*ones(ncells,p,n);
% 	par.k3 = 3*ones(p,n);
% 	par.k4 = 4*ones(p,n);
% 	par.KM = 5*ones(p,n);
% 	Pin = 7*ones(p,n);
% 	x = ones(1,(2*p*n+1)*ncells + p*n);
% 	h = 2; f1 = 1/3; f2 = 100;
% 	metab_threshold = 999*ones(1,p);

	M = zeros((p*n+1)*ncells, (p*n+1)*ncells);
	V = zeros((p*n+1)*ncells, 1);
	first_S = 1:(n*p+1):ncells*(n*p+1);
	%M_end = first_S(end) + p*n;
	first_V = 1:(n*p+1):ncells*(n*p+1);
	
	%% M
	
	count = 0;
	for ncell = 1:ncells
		%% M
		M(first_S(ncell),first_V(ncell)) = Sin;
		
		% Auxotrophy direct import
		for i = 1:n_ext		% Which of the final metabolites in the parallel cascades are externally being fed? (3 => type 1)
			M(first_S(ncell)+(n-1)*p+i, first_V(ncell)) = par.feed_rate;
		end	
		
		% Production terms in M
		len1 = [first_S(ncell)+(1:n*p)];
		len2 = [first_V(ncell)+(1:n*p)];
		M(len1,len2) = M(len1,len2) + diag(ET(ncell,:).*par.k3(:)');	
		
		% Consumption terms in M
		M(first_S(ncell),first_V(ncell)+(1:p)) = ET(ncell,:,1).*-par.k3(:,1)';
		if n > 1
			len1 = [first_S(ncell)+(1:p*(n-1))];
			len2 = [first_V(ncell)+p+(1:p*(n-1))];
			M(len1,len2) = M(len1,len2) + diag(reshape(ET(ncell,:,2:n),1,[]).*reshape(-par.k3(:,2:n),1,[]));
		end
		
		%% V
		V(first_V(ncell)) = 1;
		
		% Common Substrate step
		V(first_V(ncell)+[1:p]) = x(first_S(ncell))./( x(first_S(ncell)) + par.KM(:,1) );
		
		% Subsequent Steps in cascade
		if n > 1
			V(first_V(ncell)+p+[1:p*(n-1)]) = x(first_S(ncell)+[1:p*(n-1)])./( x(first_S(ncell)+[1:p*(n-1)]) + reshape(par.KM(:,2:n),[],1) );
		end
	end
% 	V' 
% 	M
	dx = M*V;
end