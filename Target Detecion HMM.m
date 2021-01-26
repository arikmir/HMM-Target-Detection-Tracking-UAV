% set up loop here, to iterate and predict at certain frames
  a = dir(strcat('C:\Users\Arik\Desktop\EGH400-2\DataSet\','*.png'));
 out = numel(a);
 for i = 1:out
     pngFilename = strcat('D:\Engineering\EGH400-1\hmm_viterbi\image_stream\', num2str(i), 'Result', '.png');
     image_data = imread(pngFilename);
     intensity_filter_val = 70;
     image_data(image_data < intensity_filter_val) = 0;
     image_data(image_data >= 1) = 1;
     im_array(:, :, i) = image_data;  
 end
im_array =double(im_array);
[v,h,im_array_length] = size(im_array);
frame_v = v;
frame_h = h;
%% Preprocess cropping
current_apparent_state = zeros(1, im_array_length);
j_arry = 0;
k_arry = 0;
for i = 1:im_array_length
    [j, k] = find(im_array(:,:,i) == 1);
    av_j = mean(j);
    [minValue, closestIndex] = min(abs(j - av_j.'));
    avg = floor(mean(j(closestIndex:closestIndex)));
    j = avg;
    j_arry(i) = j;
    
    av_k = mean(k);
    [minValue, closestIndex] = min(abs(k - av_k.'));
    avg = floor(mean(k(closestIndex:closestIndex)));
    k = avg;
    k_arry(i) = k;
    
    if isnan(j)
        j = [];
    end
    
    if isnan(k)
        k = [];
    end 
    image_data(image_data < intensity_filter_val) = 0;
    op1 = isempty(j);
    op2 = isempty(k);
    if k == 1 & op1 == 0 & op2 == 0
        current_apparent_state(i) =  j(1);
    elseif op1 == 0 & op2 == 0
        current_apparent_state(i) =  ((j(1)* v)-v) + k(1);
    end
end
% ((k(1) - 1)*v) + j(1)

j = 0;
k = 0;
j2 = 0;
k2 = 0;
current_apparent_state_small = zeros(1,3);
flag1 = 1;

%%%%%%%%%%%%%%%%%%%%%%% padding specs
v_s = 9;
middle_state = 41;
middle_co = 5;
range = 4;
%%%%%%%%%%%%%%%%%%%%%%%%
%% pi matrix
%N = length(image_v);
N = v_s * v_s;
pi_mat = sparse(1:N);
pi_mat(1,:) = 1/N;



%% A matrix
self_p = 7/15;
adj_p = 1/15;
image_mat = zeros(v_s,v_s);

transition_matrix = sparse(N, N);
counter = 0;
for j = 1:size(image_mat, 1)
    for i = 1:size(image_mat, 2)
        counter = counter + 1;
        image_mat_copy = image_mat;
        image_mat_copy(j, i) = self_p;    %first corner
        if j == 1 && i == 1
           image_mat_copy(j + 1, i) = adj_p;
           image_mat_copy(j, i + 1) = adj_p;
           image_mat_copy(j + 1, i + 1) = adj_p;              
        end
        if i == length(image_mat(1,:)) && j == 1 %second corner
           image_mat_copy(j, i -1) = adj_p;
           image_mat_copy(j + 1, i) = adj_p;
           image_mat_copy(j + 1, i - 1) = adj_p;              
        end
        if i == 1 && j == length(image_mat(:,1)) %third corner
           image_mat_copy(j - 1, i) = adj_p;
           image_mat_copy(j - 1, i + 1) = adj_p;
           image_mat_copy(j, i + 1) = adj_p;              
        end
        if i == length(image_mat(1,:)) && j == length(image_mat(:,1)) %fourth corner
           image_mat_copy(j, i - 1) = adj_p;
           image_mat_copy(j - 1, i - 1) = adj_p;
           image_mat_copy(j - 1, i) = adj_p;              
        end
        if j == 1 && i ~= 1 && i ~= length(image_mat(1,:))
           image_mat_copy(j, i - 1) = adj_p;
           image_mat_copy(j + 1, i - 1) = adj_p;
           image_mat_copy(j + 1, i) = adj_p;
           image_mat_copy(j + 1, i + 1) = adj_p;
           image_mat_copy(j, i + 1) = adj_p; 
        end
        if j == length(image_mat(:,1)) && i ~= 1 && i ~= length(image_mat(1,:))
           image_mat_copy(j, i - 1) = adj_p;
           image_mat_copy(j - 1, i - 1) = adj_p;
           image_mat_copy(j - 1, i) = adj_p;
           image_mat_copy(j - 1, i + 1) = adj_p;
           image_mat_copy(j, i + 1) = adj_p; 
        end
        if i == 1 && j ~= 1 && j ~= length(image_mat(:,1))
           image_mat_copy(j - 1, i) = adj_p;
           image_mat_copy(j - 1, i + 1) = adj_p;
           image_mat_copy(j, i + 1) = adj_p;
           image_mat_copy(j + 1, i + 1) = adj_p;
           image_mat_copy(j + 1, i) = adj_p; 
        end
        if i == length(image_mat(1,:)) && j ~= 1 && j ~= length(image_mat(:,1))
           image_mat_copy(j - 1, i) = adj_p;
           image_mat_copy(j - 1, i - 1) = adj_p;
           image_mat_copy(j, i - 1) = adj_p;
           image_mat_copy(j + 1, i - 1) = adj_p;
           image_mat_copy(j + 1, i) = adj_p; 
        end
        if j ~= 1 && i ~=1 && j ~= length(image_mat(:,1)) && i ~= length(image_mat(1,:))
           image_mat_copy(j - 1, i - 1) = adj_p;
           image_mat_copy(j, i - 1) = adj_p;
           image_mat_copy(j + 1, i - 1) = adj_p;
           image_mat_copy(j + 1, i) = adj_p;
           image_mat_copy(j + 1, i + 1) = adj_p; 
           image_mat_copy(j, i + 1) = adj_p; 
           image_mat_copy(j - 1, i + 1) = adj_p; 
           image_mat_copy(j - 1, i) = adj_p; 
        end
        image_mat_copy = image_mat_copy';
        image_mat_copy = image_mat_copy(:); 
        transition_matrix(:, counter) = image_mat_copy;
        image_mat_copy = image_mat;
    end
end


%% plots
figure(1)
spy(transition_matrix);



%% loop
j_n = [];
k_n = [];
estimated_frame = zeros(frame_v, frame_h);
[a, b, c] = size(im_array);
csv_output_vec = zeros(c, 1);
for i = 3:3:c % the loop
    if i ~= 1 & current_apparent_state(i) ~= 0
        current_apparent_state_small(3) = middle_state; 
        j_diff = j_arry(i) - j_arry(i-1);
        k_diff = k_arry(i) - k_arry(i-1);
        j_diff2 = j_arry(i) - j_arry(i-2);
        k_diff2 = k_arry(i) - k_arry(i-2);
        if (j_diff <= range && j_diff >= -range) && (k_diff <= range && k_diff >= -range)
            j = middle_co + j_diff;
            k = middle_co + k_diff;
            j2 = middle_co + j_diff2;
            k2 = middle_co + k_diff2;
        if k == 1
            current_apparent_state_small(2) = j;
        else
            current_apparent_state_small(2) = (k*v_s) + j;
        end
        
        if k2 == 1
            current_apparent_state_small(1) = j2;
        else
            current_apparent_state_small(1) = (k2*v_s) + j2;
        end
        end
    end 
    
im_array_new3(:,:,1) = zeros(v_s,v_s);
im_array_new3(:,:,2) = zeros(v_s,v_s);
im_array_new3(:,:,3) = zeros(v_s,v_s);
im_array_new3(5,5,3) = 1;
if j ~= 0 && k ~= 0 && j2 ~= 0 && k2 ~=0
im_array_new3(j,k,2) = 1;
im_array_new3(j2,k2,1) = 1;
end 
%%%
observation_mat = zeros(N, N);
[v,h,im_array_length] = size(im_array_new3);
for x = 1:im_array_length
    im_array_new = im_array_new3(:,:,x);
    im_array_new = im_array_new(:)';
    im_array_new_2(x,:) = im_array_new;
end
format long
present = 0;
absent = 0;
for k = 1:N
    for j = 1: im_array_length
        if im_array_new_2(j, k) == 1
            present = present+1;
        end
        if im_array_new_2(j,k) == 0
            absent = absent +1;
        end     
    end
prob_present = present/im_array_length;
prob_absent = absent/im_array_length;
obs_prob = prob_present/prob_absent;
obs_prob = double(obs_prob);
observation_mat(k,k) = obs_prob;
absent = 0;
present = 0;
obs_prob = 0;
prob_present = 0;
prob_absent = 0;
end

for c = 1:N
    for j = 1:N
       if c == j && (observation_mat(c,j)==0)
           observation_mat(c,j) = 0.0000001;
       end 
    end
end


observation_mat(observation_mat == 0) = 0.0000001;




 
if current_apparent_state_small(1) == 0 && current_apparent_state_small(2)== 0
 seq = [1:current_apparent_state_small(3) current_apparent_state_small(3) current_apparent_state_small(3)];
elseif current_apparent_state_small(1) == 0
    seq = [1:current_apparent_state_small(2) current_apparent_state_small(3) current_apparent_state_small(3) current_apparent_state_small(3)];
elseif current_apparent_state_small(2) == 0
    seq = [1:current_apparent_state_small(1) current_apparent_state_small(3) current_apparent_state_small(3) current_apparent_state_small(3)];
elseif current_apparent_state_small(1) ~= 0 && current_apparent_state_small(2) ~= 0 && current_apparent_state_small(3)~= 0
    seq = [1:current_apparent_state_small(1) current_apparent_state_small(2) current_apparent_state_small(3) current_apparent_state_small(3) current_apparent_state_small(3)];
end
%rest of the stuffs
%  [seq,states] = hmmgenerate(N,transition_matrix,observation_mat);
estimatedStates = hmmviterbi(seq, transition_matrix,observation_mat);
current_actual_state = current_apparent_state(3);%change with loop, change to i
x_arry = zeros(1,3);
y_arry = zeros(1,3);
new_j = zeros(1,3);
new_k = zeros(1,3);
x_diff = zeros(1,3);
y_diff = zeros(1,3);
for b = current_apparent_state_small(3):length(seq)%change with loop
    if estimatedStates(i) ~= 0
        x_arry(((b+1) - middle_state)) = ceil((estimatedStates(b)/v_s));
        y_arry(((b+1) - middle_state)) =  rem((estimatedStates(b)), (((x_arry(((b+1) - middle_state)))*v_s)-v_s));   
    else
        x_arry(((b+1) - middle_state)) = 0;
        y_arry(((b+1) - middle_state)) = 0;
    end
end

for k = 1:3 
    x_diff(k) = middle_co - x_arry(k);
    y_diff(k) = middle_co - y_arry(k);
end

new_j(1) = j_arry(i) + x_diff(1); %change 3 to i from main for loop
new_j(2) = j_arry(i) + x_diff(2);
new_j(3) = j_arry(i) + x_diff(3);
new_k(1) = k_arry(i) + y_diff(1); %change 3 to i from main for loop
new_k(2) = k_arry(i) + y_diff(2);
new_k(3) = k_arry(i) + y_diff(3);

final_estimated_states = zeros(1,3);
hold on 
axis([0 2050 0 2448])
ax = gca;
ax.YDir = 'reverse';
grid on
for k = 1:3
    final_estimated_states(k) = (((new_j(k) * frame_v) - frame_v) + new_k(k));
    estimated_frame(new_j(k), new_k(k)) = 1;
    scatter(new_k(k), new_j(k),'ok');
end
j_n = [j_n new_j];
k_n = [k_n new_k];
plot(k_n, j_n,'--r');
csv_output_vec(i) = final_estimated_states(1);
csv_output_vec(i-1) = final_estimated_states(2);
csv_output_vec(i-2) = final_estimated_states(3);
end



%% outside loop
writematrix(csv_output_vec, 'EstimatedStates.csv');





