trial_num = 1;  % random trial index
N=100;           % number of points
dist_max=1;     % largest pairwise distance between two different points
dist_min=0.0001;  % smallest pairwise distance between two different points
n_std = 0.00001;  % true (oracle) noise standard deviation
min_space_unit=0.00001;   % quantization step 

%%
% extraction the point locations from the solution of the proposed approach

obj_mat=[];
emd_mat=[];
p_mat=[];
for (o=[1 3 5 7 9]) % different options

    res=dlmread(strcat('./output/x_beltway_', num2str(N), '_', num2str(dist_max), '_', num2str(dist_min), '_', num2str(trial_num), '_', num2str(n_std), '_', num2str(o)));

    [res_sort, sort_idx] = sort(res, 'descend');
    sort_idx = sort_idx-1;  % location starts at 0
    loop_len = length(res)*min_space_unit;

    x_seq = sort_idx;
    x_est=[];
    for (rr=1:length(x_seq))
        if (res_sort(rr)>1e-6)                       
            x_est_tmp=[res_sort(rr) x_seq(rr)];
            x_est=[x_est; x_est_tmp];
        end
    end

    x_est(:,2)=x_est(:,2)*min_space_unit;
    [x_est_res_sort, x_est_sort_idx] = sort(x_est(:,2), 'ascend');
    x_est=x_est(x_est_sort_idx,:);

    idx = 1:size(x_est,1);
                
    while(1)
        
        x_est_choice = x_est((x_est(:,1)>0)&(x_est(:,1)<1),:);
        idx_choice = idx((x_est(:,1)>0)&(x_est(:,1)<1));
        
        merge_distance_min = max(x_est(:,2));
        merge_ri = 0;
        merge_rj = 0;
        for (ri = 1:(size(x_est_choice, 1)-1))
            for (rj = ri+1 : size(x_est_choice, 1))
                dist_loop_min = min( abs(x_est_choice(ri,2)-x_est_choice(rj,2)), loop_len - abs(x_est_choice(ri,2)-x_est_choice(rj,2)) );
                if (dist_loop_min<dist_min) 
                if (dist_loop_min<merge_distance_min) 
                    merge_distance_min = dist_loop_min;
                    merge_ri = ri;
                    merge_rj = rj;
                end
                end
            end
        end
        
        if ( (merge_ri==0) )
            break;
        else                   
            % merge_ri is the smaller one
            dist_loop_min = min( abs(x_est_choice(merge_ri,2)-x_est_choice(merge_rj,2)), loop_len - abs(x_est_choice(merge_ri,2)-x_est_choice(merge_rj,2)) );

            if (dist_loop_min==abs(x_est_choice(merge_ri,2)-x_est_choice(merge_rj,2)))
                x_est(idx_choice(merge_ri), 2) = (x_est(idx_choice(merge_ri), 1) * x_est(idx_choice(merge_ri), 2) + x_est(idx_choice(merge_rj), 1) * x_est(idx_choice(merge_rj), 2)) / (x_est(idx_choice(merge_ri), 1) + x_est(idx_choice(merge_rj), 1));
            else
                if ( x_est(idx_choice(merge_ri), 2) < x_est(idx_choice(merge_rj), 2) )
                x_est(idx_choice(merge_ri), 2) = (x_est(idx_choice(merge_ri), 1) * x_est(idx_choice(merge_ri), 2) + x_est(idx_choice(merge_rj), 1) * (x_est(idx_choice(merge_rj), 2)-loop_len) ) / (x_est(idx_choice(merge_ri), 1) + x_est(idx_choice(merge_rj), 1));
                else
                x_est(idx_choice(merge_ri), 2) = ( (x_est(idx_choice(merge_ri), 1) - loop_len) * x_est(idx_choice(merge_ri), 2) + x_est(idx_choice(merge_rj), 1) * x_est(idx_choice(merge_rj), 2) ) / (x_est(idx_choice(merge_ri), 1) + x_est(idx_choice(merge_rj), 1));
                end
            end
            x_est(idx_choice(merge_ri), 1) = x_est(idx_choice(merge_ri), 1) + x_est(idx_choice(merge_rj), 1);
            x_est(idx_choice(merge_rj), 1) = 0;
        end
                    
    end
    
    [res_sort2, sort_idx2] = sort(x_est(:,1), 'descend');
    x_est_new = x_est(sort_idx2(1:N), 2);

    
    % compute the earth mover distance
    dist_seq = [];
    for (ri=1:(length(x_est_new)-1))
        for (rj=(ri+1) : length(x_est_new))
            dist_seq = [dist_seq; norm(x_est_new(ri)-x_est_new(rj))];
            dist_seq = [dist_seq; loop_len - norm(x_est_new(ri)-x_est_new(rj))];
        end
    end
    dist_seq = [dist_seq ones(length(dist_seq), 1)];

    dist_seq_true = dlmread(strcat('./data/d_beltway_', num2str(N), '_', num2str(dist_max), '_', num2str(dist_min), '_', num2str(trial_num), '_', num2str(n_std) ));
    dist_seq_true = dist_seq_true(dist_seq_true(:,1)>0, :);
    
    dist_unique = sort(unique([dist_seq(:,1); dist_seq_true(:,1)]));
    
    dist_seq_proc = [dist_unique zeros(length(dist_unique), 1)];
    for (ri=1:size(dist_seq_proc,1))
        dist_seq_tmp = dist_seq(dist_seq(:,1)==dist_unique(ri),:);
        if (size(dist_seq_tmp,1)>0)
            dist_seq_proc(ri,2) = sum(dist_seq_tmp(:,2));
        end
    end
    
    dist_seq_true_proc = [dist_unique zeros(length(dist_unique), 1)];
    for (ri=1:size(dist_seq_true_proc,1))
        dist_seq_true_tmp = dist_seq_true(dist_seq_true(:,1)==dist_unique(ri),:);
        if (size(dist_seq_true_tmp,1)>0)
            dist_seq_true_proc(ri,2) = sum(dist_seq_true_tmp(:,2));
        end
    end
    
    emd_seq = dist_seq_proc(:,2) - dist_seq_true_proc(:,2);
    emd_seq = abs(cumsum(emd_seq));
    
    emd_seq_dist = [];
    for (ri = 2:length(dist_unique))
        emd_seq_dist = [emd_seq_dist; dist_unique(ri)-dist_unique(ri-1)];
    end
    emd_seq_dist = [emd_seq_dist; 0];
    
    emd_seq = sum(emd_seq.*emd_seq_dist);

    emd_mat = [emd_mat emd_seq];   

    % compute assignment error using Hungarian algorithm

    cost_rot = 0;
    count_rot = 0;
    p_mat_rot = [];
    for (rot_idx = 1:N)
        % rotate clock-wise through every extracted sample

        x_est = x_est_new';
        x_est = x_est-x_est(rot_idx);
        x_est(x_est<0) = x_est(x_est<0) + loop_len;
        x_ori=dlmread(strcat('./data/p_beltway_', num2str(N), '_', num2str(dist_max), '_', num2str(dist_min), '_', num2str(trial_num)));
        x_ori = x_ori';
    
        x_est_1 = x_est;
        x_est_2 = max(x_est)-x_est;
        

        cost_mat_1=zeros(N,N);
        cost_mat_2=zeros(N,N);
        for (row = 1:N)
        for (col = 1:N)
        cost_mat_1(row,col) = norm(x_est_1(:,row)-x_ori(:,col));
        cost_mat_2(row,col) = norm(x_est_2(:,row)-x_ori(:,col));
        end
        end

        [assign_1, cost_1] = munkres(cost_mat_1);
        [assign_2, cost_2] = munkres(cost_mat_2);

        cost_seq = [cost_1 cost_2];
        x_est_seq = [x_est_1; x_est_2];
        assign_seq = [assign_1' assign_2' ];
        [cost_min, idx_cost_min] = min(cost_seq);
        assign_min = assign_seq(:, idx_cost_min);
        x_est_min = x_est_seq(idx_cost_min,:);
        
        count_tmp = 0;
        for (rr=1:N)
            if (norm( x_est_min(:,rr) - x_ori(:,assign_min(rr)) )<dist_min) 
                count_tmp=count_tmp+1;
            end
        end                
        
        if (rot_idx ==1 )
            cost_rot = cost_min;
            count_rot = count_tmp;
            if (cost_1<cost_2)
                p_mat_rot = x_est_1;
            else
                p_mat_rot = x_est_2;
            end
        else
            if (cost_rot>cost_min)
                cost_rot = cost_min;
                count_rot = count_tmp;
                
                if (cost_1<cost_2)
                    p_mat_rot = x_est_1;
                else
                    p_mat_rot = x_est_2;
                end
            end
        end
    end
    obj_mat = [obj_mat count_rot];
    p_mat = [p_mat; p_mat_rot];

end


[a,b] = min(emd_mat);
rec_num = obj_mat(b);
fprintf('No. of correct points via the proposed approach: %d\n', rec_num)

p_out = p_mat(b,:);
dlmwrite(strcat('./output/p_out_beltway_', num2str(N), '_', num2str(dist_max), '_', num2str(dist_min), '_', num2str(trial_num)), p_out, 'delimiter', ' ', 'precision', 20)

figure;
scatter(x_ori, repmat(1, size(x_ori)))
hold on;
scatter(p_out, repmat(1, size(p_out)), 'r+')
legend('tru loc', 'est loc')
title(strcat('Proposed beltway reconstruction. N=', num2str(N)))

