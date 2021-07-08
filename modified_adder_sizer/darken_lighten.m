function val = darken_lighten(rgb_in, darken, intensity)	% intensity 0-1
	val = [];
	if darken == 1		% DARKEN
		val = rgb_in*(1-intensity);
	elseif darken == 0	% LIGHTEN
		val = 1 - rgb_in;
		val = rgb_in + val*intensity;
	end
end