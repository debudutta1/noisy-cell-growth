% Calculate growth rate from the inter-div_durs_compiled, by fit to an exponential growth curve
function [grow_r, grow_std] = exp_grow_rate(reps, gens, div_durs_o)
	for k = 1:reps
		div_durs = div_durs_o(randperm(length(div_durs_o)));
		div_exp_t(1) = 0;
		for i1 = 1:(gens-1)
			for i2 = 1:(2^i1)
				index = 2^i1 - 1 + i2;
				if rem(index,2) == 0	% even
					parent = index/2;	
				else
					parent = (index - 1)/2;
				end
				div_exp_t(index) = div_exp_t(parent) + div_durs(index);
			end
		end
		[j1, j2] = sort(div_exp_t/3600);
		%figure; plot(j1, log(1:length(j1)),'o');
        %% Growth curve appears to saturate at the end, because no more data points, hence the slow cells appear as saturating.
		%% Need to neglect the entire last generation of cell divisions, since there is no next generation whose fast cell divisions can continue the exponential cell division
		len_choose = fix(length(j1)*0.45);
		p = polyfit(j1(1:len_choose),log(1:len_choose),1);
		%p = polyfit(j1(1:(length(j1)-1)/2),log(1:(length(div_exp_t)-1)/2)',1);
		growth(k) = p(1);
	end
	grow_r = mean(growth);
	grow_std = std(growth);
end