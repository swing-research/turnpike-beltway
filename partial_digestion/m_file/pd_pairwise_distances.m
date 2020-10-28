C=fileread('./data/GCF_000005845.2_ASM584v2_genomic_process.fna_oneline');
C_len=length(C);



% SmaI
enzyme_idx = 1;
seq_forward = 'CCCGGG';
seq_backward = 'GGGCCC';
cutoff_pos_forward = 2;
cutoff_pos_backward = length(seq_forward) - cutoff_pos_forward - 2;

% check forward sequence locations
seq = seq_forward;
cutoff_pos = cutoff_pos_forward;

seq = char(seq);
seq_len = length(seq);
seq_str = string(seq);

pos=[0];
for (i=1:C_len-seq_len+1)
    str_C = string(C(i:i+seq_len-1));
    check = strcmp(str_C, seq_str);
    if (check~=0)
        pos = [pos i+cutoff_pos];
    end
end

% check backward sequence locations
seq = seq_backward;
cutoff_pos = cutoff_pos_backward;

seq = char(seq);
seq_len = length(seq);
seq_str = string(seq);

for (i=1:C_len-seq_len+1)
    %str_C = convertCharsToStrings(C(i:i+seq_len-1));
    str_C = string(C(i:i+seq_len-1));
    check = strcmp(str_C, seq_str);
    if (check~=0)
        pos = [pos i+cutoff_pos];
    end
end

pos = [pos C_len];
pos = sort(pos);

dlmwrite(strcat('./data/enzyme_pos_', num2str(enzyme_idx), '.mat'), pos, 'delimiter', ' ', 'precision', 10)

dist_val = [];
for (i=1:length(pos)-1)
    for (j=i+1:length(pos))
        dist_val = [dist_val pos(j)-pos(i)];
    end
end
dist_val_uq = unique(dist_val);
dist_val_uq = sort(dist_val_uq);
dist_count = [];
for (i=1:length(dist_val_uq))
    dist_count = [dist_count; dist_val_uq(i) length(dist_val(dist_val==dist_val_uq(i)))];
end

dist_count = [0 length(pos); dist_count];
dlmwrite(strcat('./data/enzyme_pw_distance_', num2str(enzyme_idx)), dist_count, 'delimiter', ' ', 'precision', 10) 




% HindIII
enzyme_idx = 2;
seq_forward = 'GGATCC';
seq_backward = 'CCTAGG';
cutoff_pos_forward = 0;
cutoff_pos_backward = length(seq_forward) - cutoff_pos_forward - 2;


% check forward sequence locations
seq = seq_forward;
cutoff_pos = cutoff_pos_forward;

seq = char(seq);
seq_len = length(seq);
seq_str = string(seq);

pos=[0];
for (i=1:C_len-seq_len+1)
    str_C = string(C(i:i+seq_len-1));
    check = strcmp(str_C, seq_str);
    if (check~=0)
        pos = [pos i+cutoff_pos];
    end
end

% check backward sequence locations
seq = seq_backward;
cutoff_pos = cutoff_pos_backward;

seq = char(seq);
seq_len = length(seq);
seq_str = string(seq);

for (i=1:C_len-seq_len+1)
    str_C = string(C(i:i+seq_len-1));
    check = strcmp(str_C, seq_str);
    if (check~=0)
        pos = [pos i+cutoff_pos];
    end
end

pos = [pos C_len];
pos = sort(pos);

dlmwrite(strcat('./data/enzyme_pos_', num2str(enzyme_idx), '.mat'), pos, 'delimiter', ' ', 'precision', 10)


dist_val = [];
for (i=1:length(pos)-1)
    for (j=i+1:length(pos))
        dist_val = [dist_val pos(j)-pos(i)];
    end
end
dist_val_uq = unique(dist_val);
dist_val_uq = sort(dist_val_uq);
dist_count = [];
for (i=1:length(dist_val_uq))
    dist_count = [dist_count; dist_val_uq(i) length(dist_val(dist_val==dist_val_uq(i)))];
end

dist_count = [0 length(pos); dist_count];
dlmwrite(strcat('./data/enzyme_pw_distance_', num2str(enzyme_idx)), dist_count, 'delimiter', ' ', 'precision', 10) 
