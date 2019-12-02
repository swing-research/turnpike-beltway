trial_num = 1;  % random trial number
rng(trial_num);

N=100;   % the number of points
dist_max=1;    % the maximum pairwise distance
dist_min = 1e-4;
L = 1.0001;    % the length of the loop, here it is dist_max + dist_min
n_std = [0 1e-5 3e-5 5e-5 7e-5 9e-5]; % the noise standard deviations under differnt noise levels


% Randomly draw N points in on a loop with length L
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

dlmwrite(strcat('./data/p_beltway_', num2str(N), '_', num2str(dist_max), '_', num2str(dist_min), '_', num2str(trial_num)), p_seq, 'delimiter', ' ', 'precision', 20)

% Compute the pairwise distances between two different points
d_seq_ori=[];
for (i=1:(length(p_seq)-1))
    for (j=(i+1):length(p_seq))
        d_seq_ori=[d_seq_ori; abs(p_seq(i)-p_seq(j))];
        d_seq_ori=[d_seq_ori; L-abs(p_seq(i)-p_seq(j))];
    end
end
d_seq_ori= sort(d_seq_ori);

% add noise to the pairwise distances
for (n=n_std)    

    d_seq = d_seq_ori + normrnd(0,n,[length(d_seq_ori) 1]);
    d_seq = abs(d_seq); % make sure the distances are positive

    d_seq=[repmat(0,N,1); d_seq];   % no need to add noise to the distance from a point to it self

    d_seq=sort(d_seq);

    d_uq_seq = unique(d_seq);
    d_out=[];
    for (i=1:length(d_uq_seq))
        d_out=[d_out; d_uq_seq(i) length(d_seq(d_seq==d_uq_seq(i)))];
    end

    dlmwrite(strcat('./data/d_beltway_', num2str(N), '_', num2str(dist_max), '_', num2str(dist_min), '_', num2str(trial_num), '_', num2str(n)), d_out, 'delimiter', ' ', 'precision', 20)
end







