function p = g_cdf(mu, sigma, x)

    % gaussian cdf
    p = 0.5*(1+erf((x-mu)/sigma/sqrt(2)));
end
