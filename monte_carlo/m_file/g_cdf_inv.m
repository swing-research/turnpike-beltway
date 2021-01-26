function x = g_cdf_inv(mu, sigma, p)

    % inverse gaussian cdf
    x = mu + sqrt(2)*sigma*erfinv(2*p-1);
end
