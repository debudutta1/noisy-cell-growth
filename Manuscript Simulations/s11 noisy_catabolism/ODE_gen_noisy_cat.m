% Function to create a Matrix form ODE for simulation of growth of cells in terms of metabolite production only. NO CROSS-FEEDING
function dx = ODE_gen_noisy_cat(t,x, ncells, p, n, Sin, par, metab_threshold, cell_dat, catabol)
	% MATLAB supposedly passes by reference by default if the variable is never modified in the function!
	ET = zeros(ncells,p,n);
	cat_E = zeros(ncells,1);
	% Interpolate the value of Total Enzyme at given time
	for cell = 1:ncells
		for parcas = 1:p
			for sercas = 1:n
				%ET(cell,parcas,sercas) = interpl_lastval(cell_dat{cell}.prot_dat{parcas,sercas}.tot_prot_profile.t,cell_dat{cell}.prot_dat{parcas,sercas}.tot_prot_profile.v,t);
				ET(cell,parcas,sercas) = interpl_lastval(cell_dat{cell}.prot_dat{parcas,sercas}.cum_prot_profile.t,cell_dat{cell}.prot_dat{parcas,sercas}.cum_prot_profile.v,t);
			end
		end
		cat_E(cell) = interpl_lastval(catabol.prot_dat{1,1}.cum_prot_profile.t,catabol.prot_dat{1,1}.cum_prot_profile.v,t);
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
	for cell = 1:ncells
		%% M
		%M(first_S(cell),first_V(cell)) = Sin;
		
		% Production terms in M
		len1 = [first_S(cell)+(1:n*p)];
		len2 = [first_V(cell)+(1:n*p)];
		M(len1,len2) = M(len1,len2) + diag(ET(cell,:).*par.k3(:)');	
		
		% Consumption terms in M
		M(first_S(cell),first_V(cell)+(1:p)) = ET(cell,:,1).*-par.k3(:,1)';
		if n > 1
			len1 = [first_S(cell)+(1:p*(n-1))];
			len2 = [first_V(cell)+p+(1:p*(n-1))];
			M(len1,len2) = M(len1,len2) + diag(reshape(ET(cell,:,2:n),1,[]).*reshape(-par.k3(:,2:n),1,[]));
		end
		
		%% V
		V(first_V(cell)) = 1;
		
		% Common Substrate step
		V(first_V(cell)+[1:p]) = x(first_S(cell)+1)./( x(first_S(cell)+1) + par.KM(:,1) );
		
		% Subsequent Steps in cascade
		if n > 1
			V(first_V(cell)+p+[1:p*(n-1)]) = x(first_S(cell)+1+[1:p*(n-1)])./( x(first_S(cell)+1+[1:p*(n-1)]) + reshape(par.KM(:,2:n),[],1) );
		end
	end
%	M0 = zeros(1, (p*n+1)*ncells);
	catab = p*par.k3(1)*cat_E(1,1) * x(first_S(1))./( x(first_S(1)) + par.KM(1) );
%	catab = p*par.k3(1)*cat_E(1,1);
	M0 = Sin - catab;
	
% 	V' 
% 	M
	M(1,1) = catab;
	dx = M*V;
	dx = vertcat(M0, dx);
end