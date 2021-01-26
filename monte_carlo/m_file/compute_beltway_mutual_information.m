function foo = estimate_beltway_mutual_information(LN_num, dist_unit, total_trial_num, output)

LN = maxNumCompThreads( str2num(LN_num) );  % set the largest number of computing threads

N=100;   % the number of points
dist_max=1;     % the maximum pairwise distance
dist_min = 5e-3;     % the minimum pairwise distance
dist_unit = str2num(dist_unit)/10000;   % the quantization step size, should be smaller than dist_min
dist_min_unit = 1/10000;                
L = dist_max+dist_min;                  % the length of the beltway loop
n_std = 0; % the noise standard deviations under differnt noise levels
total_trial_num = str2num(total_trial_num); % the number of trials



mi_trial = [];

for (trial_num = (1:total_trial_num))

    fprintf('%d\n', trial_num)

    rng(trial_num)
    % Randomly draw N points in the interval [0, dist_max]
    p_seq=[0; dist_max];   % the two outpost anchor points

    for (i=1:N-2)
        cont=1;
        while(cont==1)
            p_seq_tmp = rand(1)*dist_max;

            dist_tmp=[];
            for (j=1:size(p_seq,1))
                dist_tmp=[dist_tmp abs(p_seq_tmp-p_seq(j))];
            end

            if (min(dist_tmp)>dist_min)&&(max(dist_tmp)<dist_max)
                p_seq=[p_seq; p_seq_tmp];
                cont=0;
            end
        end
    end
    p_seq = round(p_seq/dist_unit)*dist_unit;


    L = dist_max+dist_min;

    x_ori = round(p_seq*10000);     % locations of the points in the 1d discrete domain
    K = round(L*10000);

    Ke = K;
    x_ori = x_ori + ceil(3*n_std);

    % compute the distribution of x
    x_loc = zeros(Ke,1);
    % add noise to x_ori
    for (i=1:N)
        x_idx_tmp = x_ori(i)+1;
        for (j=(x_idx_tmp-ceil(3*n_std)):(x_idx_tmp+ceil(3*n_std)))
            if ((j>=1)&&(j<=Ke))
                x_loc(j) = x_loc(j) + normcdf(j+0.5, x_idx_tmp, n_std/sqrt(2))-normcdf(j-0.5, x_idx_tmp, n_std/sqrt(2));
            end
        end 
    end

    x_loc = x_loc / sum(x_loc);

    % compute the conditonal entropy H(Y|X)
    entropy_2nd = 0;
    prob_distance = zeros(Ke,1);
    for (i=1:Ke)
        if (x_loc(i)>0)
            y_prob_tmp = zeros(Ke,1); % the distributions should be Ke+1 if we include L as a distance, here let's not do     that to complicate things
            y_prob_tmp(1:Ke-i+1) = x_loc(i:end);
            if (i>1)
                y_prob_tmp(Ke-i+2:end) = y_prob_tmp(Ke-i+2:end) + x_loc(1:i-1);
            end

            y_prob_tmp_2 = zeros(Ke,1);
            y_prob_tmp_2(1:(i-1)) = x_loc((i-1):-1:1);
            y_prob_tmp_2(i:end) = x_loc(Ke:-1:i);

            y_prob_tmp = 0.5*(y_prob_tmp + y_prob_tmp_2);

            % compute the distance distribution
            prob_distance = prob_distance + x_loc(i)*y_prob_tmp;
            
            y_prob_tmp_sum = sum(y_prob_tmp);
            if (y_prob_tmp_sum>0)
                y_prob_tmp = y_prob_tmp/y_prob_tmp_sum;
                y_prob_tmp(y_prob_tmp<=0)=[];
                entropy_2nd = entropy_2nd + x_loc(i)*sum(-y_prob_tmp.*log(y_prob_tmp));
            end
        end 
    end
    prob_distance(prob_distance<=0) = [];

    % this is the entropy of H(Y)
    entropy_1st = sum(-prob_distance.*log(prob_distance));

    % compute the mutual information I(X;Y)=H(Y)-H(Y|X)
    mi_trial = [mi_trial entropy_1st-entropy_2nd];

    fprintf('%d\t%d\n', trial_num, mi_trial(length(mi_trial)))

end

dlmwrite(strcat(output, '_mi_trial'), mi_trial, 'delimiter', ' ', 'precision', 10)

end
