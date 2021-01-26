function foo = compute_turnpike_distribution_hEh(LN_num, dist_unit, num_mc, trial_num, output)

trial_num = str2num(trial_num); % random trial number
num_mc = str2num(num_mc);       % the number of monte carlo samples
                                % note that the first 1000 gibbs samples are discarded to skip the burn-in period
                                % num_mc should be greater than 1000

LN = maxNumCompThreads( str2num(LN_num) );  % set the largest number of computing threads

rng(trial_num);     % set random seed

N=100;              % the number of points
dist_max=1;         % the maximum pairwise distance
dist_min = 5e-3;    % the minimum pairwise distance
dist_unit = str2num(dist_unit)/10000;   % the quantization step size, should be smaller than dist_min
warm_up_step = 1;   % gibbs sampling skip step

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

x_ori = round(p_seq*10000);     % the point locations
K = round(dist_max*10000)+1;


d_seq=[];
for (i=1:N)
    for (j=i:N)
        d_seq=[d_seq abs(x_ori(i)-x_ori(j))];
    end 
end

d_uq=unique(d_seq);     % unique distance values
d_uq=sort(d_uq);    


% Find the candidate locations in the 1d domain that are consisitent with the distance measurements
% this reduces the problem size, and speeds up the gibbs sampling

loc = 0:K-1;
loc_mark = zeros(K,1);
for (i=1:K)
    dist_tmp_1 = abs(loc(i)-0);
    dist_tmp_2 = abs(loc(i)-(K-1));
    if ( (min(abs(d_uq-dist_tmp_1))==0) && (min(abs(d_uq-dist_tmp_2))==0) )
        loc_mark(i)=1;
    end
end

loc_mark_num = sum(loc_mark);   % the number of candidate locations

x_loc = zeros(K,1);             % the candidate locations
x_loc(x_ori+1) = 1;

x_loc_remove = x_loc;           % remove the non-candidate locations, we only have gibbs samples in the candidate locations
x_loc_remove(loc_mark==0)=[];



% calculate the matrix E

E_seq_mat = zeros(K,K);
for (i=1:K)
    %fprintf('%d\n',i)

    E_seq_tmp_1 = zeros(K,1);
    E_seq_tmp_2 = zeros(K,1);
    E_seq_tmp_1(1:K-i+1) = x_loc(i:end);
    E_seq_tmp_2(i:end) = x_loc(1:K-i+1);
    E_seq_tmp = E_seq_tmp_1+E_seq_tmp_2;

    E_seq_mat(:,i) = E_seq_tmp;
end

% perform gibbs sampling now

Rx = repmat(1,loc_mark_num,1);
Rx(x_loc_remove==1) = -1;

[U,S,V] = svd(Rx);

% Q_x is -1*identity_matrix

Lambda = U(:,2:end)'*(diag(repmat(1,loc_mark_num,1)) - U(:,1)*U(:,1)');
Sigma_y = Lambda*U(:,2:end);
Q_y = -U(:,2:end)';
q_y = zeros(loc_mark_num,1);
L=chol(Sigma_y);
L=L';
L_T_inv = inv(L');


Q_z = L*Q_y;
q_z = zeros(loc_mark_num,1);


% initialize from z_gibbs

h_s_neg = rand(N,1);
h_s_neg = -log(h_s_neg);
h_s_neg = h_s_neg/sum(h_s_neg);

h_s_pos = rand(loc_mark_num-N,1);
h_s_pos = -log(h_s_pos);
h_s_pos = h_s_pos/sum(h_s_pos);

h_s = zeros(loc_mark_num,1);
h_s(x_loc_remove>0) = h_s_neg;
h_s((x_loc_remove<=0)) = h_s_pos;


y_h_s = U(:,2:end)'*h_s;
z_gibbs = L_T_inv*y_h_s;


z_gibbs_pre = z_gibbs;

val_sim = [];
for (i=1:num_mc)

    for (j=1:warm_up_step)
        fprintf('%d th sample, %d th warm up\n', i, j)
        z_gibbs = z_gibbs_pre;

        k=1;
        while(k<=length(z_gibbs))

            d_tmp = Q_z(k,:)';

            z_gibbs_remove = z_gibbs;
            z_gibbs_remove(k) = 0;
            Q_z_z_tmp = Q_z'*z_gibbs_remove;
            n_tmp = q_z - Q_z_z_tmp;

            d_tmp_neg = d_tmp(d_tmp<0);
            d_tmp_pos = d_tmp(d_tmp>0);

            lower_bd = -Inf;
            upper_bd = Inf;
            if (length(d_tmp_neg)>0)
                lower_bd_seq = n_tmp(d_tmp<0)./d_tmp_neg;
                lower_bd = max(lower_bd_seq);
            end
            if (length(d_tmp_pos)>0)
                upper_bd_seq = n_tmp(d_tmp>0)./d_tmp_pos;
                upper_bd = min(upper_bd_seq);
            end

            % sample from truncated gaussian via inverse transform sampling
            if (lower_bd<upper_bd)
                while(1)
                    p = rand(1);
                    z_entry = g_cdf_inv(0,1, g_cdf(0,1,lower_bd)+p*(g_cdf(0,1,upper_bd)-g_cdf(0,1,lower_bd)) );
                    if ((~isnan(z_entry))&&(~isinf(z_entry)))
                        break;
                    end
                end
                z_gibbs(k) = z_entry;
                k=k+1;
            else
                fprintf('bounds: %d\t%d\n', lower_bd, upper_bd)
                fprintf('Backtracking\n')
                z_gibbs(k-1) = z_gibbs_pre(k-1);
                k=k-1;
            end

        end

        z_gibbs_pre = z_gibbs;
    end

    % transform z_gibbs to h_s_norm
    h_s_remove = U(:,2:end)*(L'*z_gibbs);

    % to be adjusted
    h_s_remove(x_loc_remove>0) = -h_s_remove(x_loc_remove>0);

    h_s_remove_norm = h_s_remove/norm(h_s_remove);

    h_s_norm = zeros(K,1);
    h_s_norm(loc_mark==1) = h_s_remove_norm;

    val_sim_tmp = 0;
    for (j=1:K)
        val_sim_tmp = val_sim_tmp + sum(h_s_norm.*E_seq_mat(:,j))^2;
    end 

    val_sim = [val_sim val_sim_tmp];
    fprintf('val_sim %d\n', val_sim_tmp)
end

histogram(val_sim(1001:end),10)

dlmwrite(strcat(output, '_val_sim'), val_sim(1001:end), 'delimiter', ' ', 'precision', 12)  % the first 1000 gibbs samples are discarded
dlmwrite(strcat(output, '_x_ori'), x_ori, 'delimiter', ' ', 'precision', 12)
save(strcat(output, '_z_gibbs'), 'z_gibbs')

end
