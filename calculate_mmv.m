function [lambdas, types]=calculate_mmv(lms1, cm1, n_points)
% function lambdas=calcmmv(LMS1, illum, cones, wavelengths, nPoints)
% finds nPoints random points metameric to lms1 under illum
% in the form of spectra described by 5 transition wavelengths
% each row of "lambdas" describes one spectrum
% type is determined by order:
% for Type I: L1<=L2<=...
% for Type II: L1>=L2>=...
% (type is then automatically dealt with by all other functions)

% display(sprintf('calcmmv: generating %d points on MMV surface...',nPoints));
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization
    maxError=0.005; % maximum error for metamers in LMS space

    % The functions to minimize
    % we want optimal n-transition reflectances metameric to lms1 under cm1
    fT1=@(WLs) sqrt( sum((cm1.ntr2lms_r(sort(WLs),1)-lms1).^2,2) );
    fT2=@(WLs) sqrt( sum((cm1.ntr2lms_r(sort(WLs),2)-lms1).^2,2) );
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find metamers by optimization
n_points=round(n_points/2);
% Type I
% disp('Finding Type 1 metamers...');

n_transitions=size(lms1,2)*2-1;
% n_transitions=7;
lambdas=zeros(n_points*2,n_transitions); % the output 5-transition spectra metameric to LMS1
lb=lambdas(1,:)*0;
ub=lb+1;
options=optimset('MaxFunEvals',50000,'MaxIter',50000);
 
fprintf('\n Finding % i points: ',n_points*2);
parfor i=1:n_points
%     fprintf('%i..',i);
    
    if mod(i,100)==0
        disp(i);
    end

    x0=rand(1,n_transitions);
    % solve for initial metamer
    xn=optimize(fT1, x0, lb, ub, [], [], [], [], [], [], options);
    xbest=xn;
    uff=0;
    % now iterate until the error is small enough
    while ~(fT1(xbest)<=maxError)
        uff=uff+1;
        x0=rand(1,n_transitions);
        xn=optimize(fT1, x0, lb, ub, [], [], [], [], [], [], options);
        if fT1(xn)<fT1(xbest)
            xbest=xn;
        end
    end
    lambdas(i,:)=sort(xbest);
end

% Type II
% disp('Finding Type 2 metamers...');

parfor i=(n_points+1):(n_points*2)
%     fprintf('%i..',i);
    if mod(i,100)==0
        disp(i);
    end

    x0=rand(1,n_transitions);
    % solve for initial metamer
    xn=optimize(fT2, x0, lb, ub, [], [], [], [], [], [], options);
    xbest=xn;
    uff=0;
    % now iterate until the error is small enough
    while ~(fT2(xbest)<=maxError)
        uff=uff+1;
        x0=rand(1,n_transitions);
        xn=optimize(fT2, x0, lb, ub, [], [], [], [], [], [], options);
        if fT2(xn)<fT2(xbest)
            xbest=xn;
        end
    end
    lambdas(i,:)=sort(xbest,'descend');
end
types=( (1:n_points*2) > n_points )+1;
fprintf('\n Done! \n');

lambdas(:)=cm1.deparametrize(lambdas(:));
elapsed=toc;
% display(sprintf('calcmmv finished, time passed: %.2f seconds \n',elapsed));
