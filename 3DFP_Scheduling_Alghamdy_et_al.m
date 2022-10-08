clear all; clc
% TIGHT DEADLINES LIMDA 1000
p_set = 1; % machine setup time
p_vt = 0.030864; % Time spent to form per unit volume of material
% volume for each part
P_vol = [826.08; 952.60; 71.91; 703.08; 272.92; 125.70; ...
    1142.25; 121.82; 315.00; 102.83; 214.79; 124.66];
p_ma = 900; % Machine area
% Area for each part
P_area = [209.06; 550.11; 23.63; 99.53; 56.85; 50.02; ...
    435.66; 84.97; 48.27; 122.62; 178.34; 134.08];
% Height for each part 
P_h = [6.90; 26.04; 15.97; 17.04; 27.94; 17.38; ...
    11.81; 2.67; 17.13; 4.27; 2.18; 6.48]; 
% Time spent for powder-layering, which is repeated for each layer based on 
% the highest part produced in the job
p_ht = 0.7;   
I = 12; % total number of parts
J = 3; % total number of jobs
% the lower bound for each job
p_LB = [0; 0; 0]; 
% the upper bound for each job
p_UB = [14; 92; 188];



% Hyper-parameter for the job utilisation constraint

P_ddline = [90; 180; 187; 190; 188; 200; 92; 96; 189; 8; 10; 98]; % deadline for each part

% the regularization parameter to control the weight of delay 
% minimization relative to PT minimization 
LAMBDA = [0, 1e-3, 1e-2, 1e-1, 0.5, 1, 10];

PT_total = zeros(1, length(LAMBDA));
DV_total = zeros(1, length(LAMBDA));
%% the automated part of the code
for itr = 1:length(LAMBDA)
lambda = LAMBDA(itr); 

%% the cost function
f = [repmat(p_vt.*P_vol', 1, J), ...
    p_set.*ones(1, J), ...
    p_ht.*ones(1, J), ...
    zeros(1, J), lambda.*ones(1, J)];
%% equality constraints Eq. 4 
Aeq = zeros(I, I*J + 4*J);
for i = 1:I
    for j = 1:J
        Aeq(i, I*(j-1) + 1 + (i - 1)) = 1;
    end
end
beq = ones(I, 1);
%% eq# 0 (define zj's)
A0 = [eye(I*J, I*J), zeros(I*J, 4*J)];
for j = 1:J
    A0(I*(j-1)+1:I*j, (I*J)+j) = -1;
end
b0 = zeros(I*J, 1);
%% eq# 5
A5 = zeros(J, I*J+ 4*J);
for j = 1:J
   A5(j, (j-1)*I+1: j*I)  = P_area';
end
b5 = [ones(J,1).*p_ma];
%% eq# 6: the job utilization equation (we decided to omit it)
% A6 = zeros((J-1), (I*J)+(4*J));
% for j = 1:J-1
%     A6(j, (j-1)*I+1: j*I) = -p_psi.*ones(1, I);
%     A6(j, j*I+1: (j+1)*I) = ones(1, I);
% end
% b6 = zeros(J-1,1);
%% eq# 9 (to deal with max in the cost function of PT)
A9 = [diag(repmat(P_h, [J, 1])), zeros(I*J, 4*J)];
for j = 1:J
    A9(I*(j-1)+1:I*j, (I*J)+J+j) = -1;
end
b9 = zeros(I*J, 1);

%% Delay time Eq. (Eq. 6)
Adt = zeros(J, I*J + 4*J);
for j = 1:J
    Adt(j, 1: j*I) = repmat(p_vt.*P_vol', 1, j);
    Adt(j, I*J+J+1:I*J+J+j) = p_set.*ones(1, j);
    Adt(j, I*J+J+1:I*J+J+j) = p_ht.*ones(1, j);
    Adt(j, I*J+2*J+j) = -1;
    Adt(j, I*J+3*J+j) = -1;
end
bdt = zeros(1*J, 1);
%% final constraints 
A = [A0; A5; A9; Adt];
b = [b0; b5; b9; bdt];
% final problem settings
intcon = 1:I*J+J;
lb = [zeros(I*J+2*J,1); p_LB; zeros(J, 1)]; 
ub = [ones(I*J + J,1); ones(J, 1)*Inf; p_UB; ones(J, 1)*Inf]; % Enforces x(3) is binary
%% optimization toolbox
tic
x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
toc
xj = [];
for j = 1:J
    xj = [xj, x((j-1)*I+1:j*I)];
end
%% formating and printing results
X_final = cell(1, J);
PT = cell(1, J);
for j = 1:J
    [partIdx, ~] = find(xj(:,j) > 0.999);
    X_final{j} = partIdx;
    PTIdx = [(j-1)*I+1:j*I, I*J+j, I*J+J+j];
    PT{j} = f(PTIdx)*x(PTIdx);
end

f = [p_vt.*P_vol', p_vt.*P_vol', p_vt.*P_vol', p_set.*ones(1, J), p_ht.*ones(1, J), ...
    zeros(1, J), lambda.*ones(1, J)];

sumPT = 0;
for j = 1:J
    disp(['***** The ID of parts in Job number ', num2str(j), ' are: *****'])
    disp(X_final{j}')
    disp(['****** The PT of Job number ', num2str(j), ' is'])
    disp(PT{j})
    sumPT = sumPT + PT{j};
end
%% presentable results 
% printing the total PT time
disp(['****** The total PT for all Jobs is']) 
disp(sumPT)
PT_total(itr) = sumPT;

% printing the deadline violation
ddlineViolation = zeros(1, I);
absPTj = 0;
for j = 1:J
    absPTj = absPTj + PT{j};
    for i = 1:length(X_final{j})
%        disp(['job id is ', num2str(j), '*** part id is ', num2str(X_final{j}(i))])
        ddlineViolation(X_final{j}(i)) = max(0, absPTj - P_ddline(X_final{j}(i)));
    end
end

disp('*** The deadline violation for each parts is ')
disp(ddlineViolation)
disp('*** The TOTAL deadline violation for all parts is ')
disp(sum(ddlineViolation))

DV_total(itr) = sum(ddlineViolation);
end
%% lambda plot
x_axis = (LAMBDA);
y = (DV_total);


xt = x_axis;
xv = 1:length(xt);
yv = interp1(x_axis, y, xt, 'linear', 'extrap');             % Extrapolate The ‘3500’ Value
figure(1)
yyaxis left
plot(xv, yv, 'b--o')                                            % Plot Interpolated Values Against ‘xv’
grid on
set(gca, 'XTickLabel',xt)                               % Label Ticks As ‘xt’
ylabel('Total Deadline Violation (hrs)')
xlabel('Lambda Value')

% figure(2)
yyaxis right
plot(xv,PT_total, 'r--*');
grid
set(gca, 'XTickLabel',xt)                               % Label Ticks As ‘xt’
ylabel('Total Process Time (hrs)')
xlabel('Lambda Value')
ylim([180 190])
grid on