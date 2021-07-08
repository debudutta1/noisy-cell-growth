% Binomial distribution approximated as truncated normal
%    If X ~ B(n, p) and if n is large and/or p is close to ½, then X is approximately N(np, npq), (where q = 1 - p).

function y = bino_normal_rnd(n)
	s = size(n);
	for i = 1:s(1)
		for j = 1:s(2)	
			if n(i,j) ~= 0
                pd(i,j) = makedist('Normal', n(i,j)/2, n(i,j)/4);
                pd(i,j) = truncate(pd(i,j),0,n(i,j));
                y(i,j) = floor(random(pd(i,j), 1));		% Use lower instead of round, since round could round up to be above the max!
            else
                y(i,j) = 0;
            end
		end
	end
end