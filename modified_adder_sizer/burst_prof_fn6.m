% Function assumes that all the enzymes are on the same operon and thus creates a single stochastic profile and assumes that the subsequent enzymes expressed the same as the first but with a time delay equal to the size of the first gene transcribed.
function y = burst_prof_fn6(clock, n, T, transcribe_time, t_on, t_off)%, extend, burst_end)	% n-cascade size, clock - start time, T - time profile length	
	% Exponential distribution D1 to derive the duration of transcription ON time
	%t_on = 6*60; % On Average transcription ON time in minutes, converted to seconds
	% Exponential distribution D2 to derive the duration of transcription OFF time
	%t_off = 37*60; % On Average transcription OFF time in minutes, converted to seconds

	start_t = [];
	stop_t = [];
	
	% Since burst durations are distributed exponentially, which is a memory less distribution the probability of burst durations remains the same, wherever we start measuring. Thus, it is not necessary to know where the last burst ended!
	time = clock;
	r = 0;
	
	% Create burst profile
	while time <= clock+T
		r = random('exp', t_off);
		time = time + r;			% Start of burst
		if time <= clock+T		% Continue the burst record only if next burst event is before clock+T
			start_t = horzcat(start_t, time);	% Add to burst period start vector list
			r = random('exp', t_on); 	% Duration of the burst
			time = time + r;			% End of burst
			stop_t = horzcat(stop_t, time);	% Add to burst period start vector list
		end	
	end	
	
	% Output variables
	y.start = start_t;
	y.stop = stop_t;
    if ~isempty(stop_t)
        y.burst_end = stop_t(end);
    else
        y.burst_end = clock+T;
    end
end