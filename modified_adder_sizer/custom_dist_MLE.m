% Fit to custom distribution using MLE

load('readmissiontimes.mat');

%Define a custom probability density and cumulative distribution function.

custpdf = @(data,lambda) lambda*exp(-lambda*data);
custcdf = @(data,lambda) 1-exp(-lambda*data);

%Estimate the parameter, lambda, of the custom distribution for the censored sample data.
phat = mle(ReadmissionTime,'pdf',custpdf,'cdf',custcdf,'start',0.05,'Censoring',Censored)


%==================

% Fréchet distribution

% 3 Parameters: k - shape (alpha), sig - scale (beta), mu - location (thetha). k > 0, sig > 0

frechet_pdf = @(x, k, sig, mu) (k/sig).*((x-mu)/sig).^(-k-1).*exp(-((x-mu)/sig).^(-k));
frechet_cdf = @(x, k, sig, mu) exp(-((x-mu)/sig).^(-k));
phat = mle(doubDat_std, 'pdf', frechet_pdf, 'cdf', frechet_cdf, 'start', [0.1, 1, -0.5]);

% 2 Parameter: k - shape (a), sig - scale (b). k > 0, sig > 0, x>0
frechet_pdf = @(x, k, sig) (k/sig)*((sig/x).^(k+1)).*exp(-(sig/x).^k);
frechet_cdf = @(x, k, sig) exp(-(sig/x).^k);

frechet_pdf = @(x, k, sig) (k/sig)*((sig/x')^(k+1))*exp(-(sig/x')^k);
frechet_cdf = @(x, k, sig) exp(-(sig/x')^k);

phat = mle(doubDat_std, 'pdf', frechet_pdf, 'cdf', frechet_cdf, 'Lowerbound',[0, 0], 'start', [0.1, 1]);
phat = mle(doubDat_std, 'pdf', frechet_pdf, 'start', [0.1, 1]);

% Log PDF

frechet_logpdf = @(x, k, sig, mu) log(k)-log(sig)+(k+1)*log(sig)-(k+1)*log(x-mu)-(sig/(x-mu)).^k;
phat = mle(doubDat_std, 'logpdf', frechet_logpdf,'start', [0.1, 1, -0.5], 'Lowerbound',[0, 0, -2]);

