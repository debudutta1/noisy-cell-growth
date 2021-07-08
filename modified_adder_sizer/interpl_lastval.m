% Function to return the last queried value for single xq query
function vq = interpl_lastval(x,v,xq)
	index = find(x > xq);
	if ~isempty(index)
		% If the queried value is large than the last value in the input, use the last value in the input	
		vq = v(index(1)-1);
	else
		% If the queried value is larger than input range, return the last value in range
		vq = v(end);
	end
end