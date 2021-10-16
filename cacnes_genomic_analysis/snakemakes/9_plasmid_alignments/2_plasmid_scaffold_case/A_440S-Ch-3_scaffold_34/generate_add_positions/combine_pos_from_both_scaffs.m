load('other_positions_3.mat')
op3=op;
load('other_positions_4.mat')
op4=op;
op=[op3;numel(op3)+op4];
isequal(op',1:1:numel(op3)+numel(op4));
save('other_positions.mat','op')