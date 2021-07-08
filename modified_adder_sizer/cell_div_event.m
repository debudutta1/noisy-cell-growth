% Events function for the ODE generated
function [position,isterminal,direction] = cell_div_event(t,y,ncells,p,n,metab_threshold, each_cell)
	value = true(ncells,1);
	position = ones(ncells,1);
    isterminal = ones(ncells,1);
    direction = zeros(ncells,1);
	for k = 1:ncells
		for m = 1:p
			value(k) = value(k) & (y(each_cell*k-p+m) > metab_threshold(m));
		end
		position(k) = single(~value(k));
    end
end