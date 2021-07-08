% Function to create a compiled simulated cell growth record!
function rec = compile_cell_rec(n,p,birth_t, div_t, cell_dat, sol2, mrna_end, prot_end, ES_ratio, parent, mrna_beg, prot_beg, detailed_Rec)		
% Input birth and division times, cell_dat for enzyme profile, and sol2 for metabolite profile, mrna_end/prot_end - interpolated values at division

	if detailed_Rec == 1
		rec = struct('mrna_prof',{{}},'prot_prof',{{}},'metab_prof',[],'birth_t',0,'div_t',0, 'ES_ratio',[]);		
		for i = 1:p
			for j = 1:n
				% Joined mrna profile
				ind = (cell_dat.mrna_dat{i,j}.cum_rna_profile.t >= birth_t & cell_dat.mrna_dat{i,j}.cum_rna_profile.t < div_t);
				rec.mrna_prof{i,j}.t = cell_dat.mrna_dat{i,j}.cum_rna_profile.t(ind);
				rec.mrna_prof{i,j}.v = cell_dat.mrna_dat{i,j}.cum_rna_profile.v(ind);
				rec.mrna_prof{i,j}.t = cat(2, rec.mrna_prof{i,j}.t, div_t);
				rec.mrna_prof{i,j}.v = cat(2, rec.mrna_prof{i,j}.v, mrna_end(i,j));
				
				% Joined prot profile
				ind = (cell_dat.prot_dat{i,j}.cum_prot_profile.t >= birth_t & cell_dat.prot_dat{i,j}.cum_prot_profile.t < div_t);
				rec.prot_prof{i,j}.t = cell_dat.prot_dat{i,j}.cum_prot_profile.t(ind);
				rec.prot_prof{i,j}.v = cell_dat.prot_dat{i,j}.cum_prot_profile.v(ind);
				rec.prot_prof{i,j}.t = cat(2, rec.prot_prof{i,j}.t, div_t);
				rec.prot_prof{i,j}.v = cat(2, rec.prot_prof{i,j}.v, prot_end(i,j));
			end
		end
		rec.metab_prof = struct('x',[],'y',[]);
		for i = length(sol2):-1:1
			if sol2{i}.x >= birth_t & sol2{i}.x <= div_t
				rec.metab_prof.x = cat(2,sol2{i}.x, rec.metab_prof.x);
				rec.metab_prof.y = cat(2,sol2{i}.y, rec.metab_prof.y);
			end
		end
	end
	rec.mrna_beg = mrna_beg;
	rec.prot_beg = prot_beg;
	rec.mrna_end = mrna_end;
	rec.prot_end = prot_end;
	rec.ES_ratio = ES_ratio;
	rec.parent = parent;
	rec.birth_t = birth_t;
	rec.div_t = div_t;
end