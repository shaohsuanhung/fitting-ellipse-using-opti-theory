%% 2DME20 PROJECT TU/E
%% SHAO HSUAN HUNG & FLORIS NABER 

%% Question 3 %%
f = imread('ClearMoon.jpg'); % load image
f_gray = rgb2gray(f);        % convert to grayscale

% smooth image
sigma_gauss = 10;
f_smooth = imgaussfilt(f_gray,sigma_gauss);

% edge detection
threshold = [0.4, 0.9]; sigma_canny = 3;
f_feature_map = edge(f_smooth,"Canny",threshold,sigma_canny);

% feature set
[row, col] = find(f_feature_map==1); %find coords of edges
sz = size(f_feature_map);
f_feature_set = [row-sz(1)/2,col-sz(2)/2]; %(y,x) y=[-.5*Y, .5*Y], x=[-.5*X, .5*X]
clear row; clear col;

p = figure('Position',[20,50,800,500]);
subplot(1,3,1)
imshow(f_gray); colormap gray
title("Grayscale image")
subplot(1,3,2)
imshow(f_smooth); colormap gray
title("Smoothed image")
subplot(1,3,3)
imshow(f_feature_map); colormap gray
title("Edges detected in image")
saveas(p,'q3.eps','epsc')

%% linear least squares residual distance
A = [f_feature_set(:,1).^2 2*f_feature_set(:,1).*f_feature_set(:,2) f_feature_set(:,2).^2 2*f_feature_set(:,1) 2*f_feature_set(:,2)];
b = ones(length(f_feature_set),1);
x = (inv(A'*A))*(A'*b); %solution of least squares
%solve system of equations
e_ls = (x(2)^2-x(1)*x(3))/(x(1)*x(3)-x(2)^2+x(1)*x(5)^2+x(3)*x(4)^2-2*x(2)*x(4)*x(5)); % solve for e
parameters = -e_ls.*x; % ellipse parameters
Q_ls = [parameters(1) parameters(2);parameters(2) parameters(3)]; d_ls = [parameters(4); parameters(5)]; %ellipse parameters Q and d
c_ls = -inv(Q_ls)*d_ls; L_ls = chol(Q_ls); E_ls = inv(L_ls); %compete c and E for plot
X = 0:0.01:2*pi;
coords = [cos(X);sin(X)];
ellipse1 = c_ls + E_ls*coords;
ellipse1 = [ellipse1(1,:)+sz(1)/2; ellipse1(2,:)+sz(2)/2];
P_tot1 = sum((A*parameters+b*e_ls).^2);

% try some other images

f2 = rgb2gray(imread('eggplant.png')); f2_2 = edge(imgaussfilt(f2,sigma_gauss),"Canny",threshold,sigma_canny);
[row, col] = find(f2_2==1); sz2 = size(f2_2); f_feature_set2 = [row-sz2(1)/2,col-sz2(2)/2]; clear row; clear col;

f3 = rgb2gray(imread('lizzard.png')); f3_2 = edge(imgaussfilt(f3,sigma_gauss),"Canny",threshold,sigma_canny);
[row, col] = find(f3_2==1); sz3 = size(f3_2); f_feature_set3 = [row-sz3(1)/2,col-sz3(2)/2]; clear row; clear col;

A = [f_feature_set2(:,1).^2 2*f_feature_set2(:,1).*f_feature_set2(:,2) f_feature_set2(:,2).^2 2*f_feature_set2(:,1) 2*f_feature_set2(:,2)];
b = ones(length(f_feature_set2),1);
x = (inv(A'*A))*(A'*b); %solution of least squares
%solve system of equations
e_ls = (x(2)^2-x(1)*x(3))/(x(1)*x(3)-x(2)^2+x(1)*x(5)^2+x(3)*x(4)^2-2*x(2)*x(4)*x(5)); % solve for e
parameters = -e_ls.*x; % ellipse parameters
Q_ls = [parameters(1) parameters(2);parameters(2) parameters(3)]; d_ls = [parameters(4); parameters(5)]; %ellipse parameters Q and d
c_ls = -inv(Q_ls)*d_ls; L_ls = chol(Q_ls); E_ls = inv(L_ls); %compete c and E for plot
ellipse2 = c_ls + E_ls*coords;
ellipse2 = [ellipse2(1,:)+sz2(1)/2; ellipse2(2,:)+sz2(2)/2];
P_tot2 = sum((A*parameters+b*e_ls).^2);

A = [f_feature_set3(:,1).^2 2*f_feature_set3(:,1).*f_feature_set3(:,2) f_feature_set3(:,2).^2 2*f_feature_set3(:,1) 3*f_feature_set3(:,2)];
b = ones(length(f_feature_set3),1);
x = (inv(A'*A))*(A'*b); %solution of least squares
%solve system of equations
e_ls = (x(2)^2-x(1)*x(3))/(x(1)*x(3)-x(2)^2+x(1)*x(5)^2+x(3)*x(4)^2-2*x(2)*x(4)*x(5)); % solve for e
parameters = -e_ls.*x; % ellipse parameters
Q_ls = [parameters(1) parameters(2);parameters(2) parameters(3)]; d_ls = [parameters(4); parameters(5)]; %ellipse parameters Q and d
c_ls = -inv(Q_ls)*d_ls; L_ls = chol(Q_ls); E_ls = inv(L_ls); %compete c and E for plot
ellipse3 = c_ls + E_ls*coords;
ellipse3 = [ellipse3(1,:)+sz3(1)/2; ellipse3(2,:)+sz3(2)/2];
P_tot3 = sum((A*parameters+b*e_ls).^2);

p = figure('Position',[20,50,1000,500]);
subplot(2,3,1)
imshow(f_feature_map); colormap gray; hold on
plot(ellipse1(2,:),ellipse1(1,:),'r',LineWidth=2); hold off
title("Ellipse fitted on edges of moon")
subplot(2,3,2)
imshow(f2_2); colormap gray; hold on
plot(ellipse2(2,:),ellipse2(1,:),'r',LineWidth=2); hold off
title("Ellipse fitted on edges of eggplant")
subplot(2,3,3)
imshow(f3_2); colormap gray; hold on
plot(ellipse3(2,:),ellipse3(1,:),'r',LineWidth=2); hold off
title("Ellipse fitted on edges of lizard")
subplot(2,3,4)
imshow(f_gray); colormap gray; hold on
plot(ellipse1(2,:),ellipse1(1,:),'r',LineWidth=2); hold off
title("Ellipse fit overlayed on moon")
subplot(2,3,5)
imshow(f2); colormap gray; hold on
plot(ellipse2(2,:),ellipse2(1,:),'r',LineWidth=2); hold off
title("Ellipse fit overlayed on eggplant")
subplot(2,3,6)
imshow(f3); colormap gray; hold on
plot(ellipse3(2,:),ellipse3(1,:),'r',LineWidth=2); hold off
title("Ellipse fit overlayed on lizard")
saveas(p,'q6.eps','epsc')


%% linear program
A = [f_feature_set(:,1).^2 2*f_feature_set(:,1).*f_feature_set(:,2) f_feature_set(:,2).^2 2*f_feature_set(:,1) 2*f_feature_set(:,2)];
b = ones(length(f_feature_set),1);
R = optimvar('R','LowerBound',0); % variable in objective function
for i=1:length(b)
    d_res_pos_n(i,:) = string((['d_res_pos', num2str(i)]));
    d_res_neg_n(i,:) = string((['d_res_neg', num2str(i)]));
end
d_res_pos = optimvar('d_res_pos',d_res_pos_n,'LowerBound',0); %create variables for d_res of all measurements
d_res_neg = optimvar('d_res_neg',d_res_neg_n,'LowerBound',0);
clear d_res_pos_n; clear d_res_neg_n;
q11_e = optimvar('q11_e','LowerBound',0); q12_e = optimvar('q12_e'); q22_e = optimvar('q22_e','LowerBound',0); %Q/-e variables
d1_e = optimvar('d1_e'); d2_e = optimvar('d2_e'); %d/-e variables 
linprob = optimproblem('Objective', R);
for i=1:length(b)
    linprob.Constraints.(['econs', num2str(i)]) = d_res_pos(i) - d_res_neg(i) == ...
        A(i,1)*q11_e + A(i,2)*q12_e + A(i,3)*q22_e + A(i,4)*d1_e + A(i,5)*d2_e - b(i); %d_res = Ax+b
    linprob.Constraints.(['cons', num2str(i)]) = d_res_pos(i) + d_res_neg(i) <= R; %abs(d_res) should be smaller than R
end
linsol = solve(linprob); %solve the program
x = [linsol.q11_e; linsol.q12_e; linsol.q22_e; linsol.d1_e; linsol.d2_e];

e_lp = (x(2)^2-x(1)*x(3))/(x(1)*x(3)-x(2)^2+x(1)*x(5)^2+x(3)*x(4)^2-2*x(2)*x(4)*x(5)); % solve for e
parameters = -e_lp.*x; % ellipse parameters
Q_lp = [parameters(1) parameters(2);parameters(2) parameters(3)]; d_lp = [parameters(4); parameters(5)]; %ellipse parameters Q and d
c_lp = -inv(Q_lp)*d_lp; L_lp = chol(Q_lp); E_lp = inv(L_lp); %compete c and E for plot
ellipse1 = c_lp + E_lp*coords;
ellipse1 = [ellipse1(1,:)+sz(1)/2; ellipse1(2,:)+sz(2)/2];
P_max1 = -e_lp*linsol.R; % P_max1 = max(abs(A*parameters+b*e_lp))

% try it for some other images
A = [f_feature_set2(:,1).^2 2*f_feature_set2(:,1).*f_feature_set2(:,2) f_feature_set2(:,2).^2 2*f_feature_set2(:,1) 2*f_feature_set2(:,2)];
b = ones(length(f_feature_set2),1);
for i=1:length(b)
    d_res_pos_n(i,:) = string((['d_res_pos', num2str(i)]));
    d_res_neg_n(i,:) = string((['d_res_neg', num2str(i)]));
end
d_res_pos = optimvar('d_res_pos',d_res_pos_n,'LowerBound',0); %create variables for d_res of all measurements
d_res_neg = optimvar('d_res_neg',d_res_neg_n,'LowerBound',0);
clear d_res_pos_n; clear d_res_neg_n;
linprob2 = optimproblem('Objective', R);
for i=1:length(b)
    linprob2.Constraints.(['econs', num2str(i)]) = d_res_pos(i) - d_res_neg(i) == ...
        A(i,1)*q11_e + A(i,2)*q12_e + A(i,3)*q22_e + A(i,4)*d1_e + A(i,5)*d2_e - b(i); %d_res = Ax+b
    linprob2.Constraints.(['cons', num2str(i)]) = d_res_pos(i) + d_res_neg(i) <= R; %abs(d_res) should be smaller than R
end
linsol = solve(linprob2); %solve the program
x = [linsol.q11_e; linsol.q12_e; linsol.q22_e; linsol.d1_e; linsol.d2_e];
e_lp = (x(2)^2-x(1)*x(3))/(x(1)*x(3)-x(2)^2+x(1)*x(5)^2+x(3)*x(4)^2-2*x(2)*x(4)*x(5)); % solve for e
parameters = -e_lp.*x; % ellipse parameters
Q_lp = [parameters(1) parameters(2);parameters(2) parameters(3)]; d_lp = [parameters(4); parameters(5)]; %ellipse parameters Q and d
c_lp = -inv(Q_lp)*d_lp; L_lp = chol(Q_lp); E_lp = inv(L_lp); %compete c and E for plot
ellipse2 = c_lp + E_lp*coords;
ellipse2 = [ellipse2(1,:)+sz2(1)/2; ellipse2(2,:)+sz2(2)/2];
P_max2 = -e_lp*linsol.R; % P_max2 = max(abs(A*parameters+b*e_lp))

A = [f_feature_set3(:,1).^2 2*f_feature_set3(:,1).*f_feature_set3(:,2) f_feature_set3(:,2).^2 2*f_feature_set3(:,1) 2*f_feature_set3(:,2)];
b = ones(length(f_feature_set3),1);
for i=1:length(b)
    d_res_pos_n(i,:) = string((['d_res_pos', num2str(i)]));
    d_res_neg_n(i,:) = string((['d_res_neg', num2str(i)]));
end
d_res_pos = optimvar('d_res_pos',d_res_pos_n,'LowerBound',0); %create variables for d_res of all measurements
d_res_neg = optimvar('d_res_neg',d_res_neg_n,'LowerBound',0);
clear d_res_pos_n; clear d_res_neg_n;
linprob3 = optimproblem('Objective', R);
for i=1:length(b)
    linprob3.Constraints.(['econs', num2str(i)]) = d_res_pos(i) - d_res_neg(i) == ...
        A(i,1)*q11_e + A(i,2)*q12_e + A(i,3)*q22_e + A(i,4)*d1_e + A(i,5)*d2_e - b(i); %d_res = Ax+b
    linprob3.Constraints.(['cons', num2str(i)]) = d_res_pos(i) + d_res_neg(i) <= R; %abs(d_res) should be smaller than R
end
linsol = solve(linprob3); %solve the program
x = [linsol.q11_e; linsol.q12_e; linsol.q22_e; linsol.d1_e; linsol.d2_e];
e_lp = (x(2)^2-x(1)*x(3))/(x(1)*x(3)-x(2)^2+x(1)*x(5)^2+x(3)*x(4)^2-2*x(2)*x(4)*x(5)); % solve for e
parameters = -e_lp.*x; % ellipse parameters
Q_lp = [parameters(1) parameters(2);parameters(2) parameters(3)]; d_lp = [parameters(4); parameters(5)]; %ellipse parameters Q and d
c_lp = -inv(Q_lp)*d_lp; L_lp = chol(Q_lp); E_lp = inv(L_lp); %compete c and E for plot
ellipse3 = c_lp + E_lp*coords;
ellipse3 = [ellipse3(1,:)+sz3(1)/2; ellipse3(2,:)+sz3(2)/2];
P_max3 = -e_lp*linsol.R; % P_max3 = max(abs(A*parameters+b*e_lp))

p = figure('Position',[20,50,1000,500]);
subplot(2,3,1)
imshow(f_feature_map); colormap gray; hold on
plot(ellipse1(2,:),ellipse1(1,:),'r',LineWidth=2); hold off
title("Ellipse fitted on edges of moon")
subplot(2,3,2)
imshow(f2_2); colormap gray; hold on
plot(ellipse2(2,:),ellipse2(1,:),'r',LineWidth=2); hold off
title("Ellipse fitted on edges of eggplant")
subplot(2,3,3)
imshow(f3_2); colormap gray; hold on
plot(ellipse3(2,:),ellipse3(1,:),'r',LineWidth=2); hold off
title("Ellipse fitted on edges of lizard")
subplot(2,3,4)
imshow(f_gray); colormap gray; hold on
plot(ellipse1(2,:),ellipse1(1,:),'r',LineWidth=2); hold off
title("Ellipse fit overlayed on moon")
subplot(2,3,5)
imshow(f2); colormap gray; hold on
plot(ellipse2(2,:),ellipse2(1,:),'r',LineWidth=2); hold off
title("Ellipse fit overlayed on eggplant")
subplot(2,3,6)
imshow(f3); colormap gray; hold on
plot(ellipse3(2,:),ellipse3(1,:),'r',LineWidth=2); hold off
title("Ellipse fit overlayed on lizard")
saveas(p,'q7.eps','epsc')