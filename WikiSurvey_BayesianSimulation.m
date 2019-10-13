function WikiSurvey_BayesianSimulation
%----------------------------------------------------------------------------------------------
% Bayesian analysis
%----------------------------------------------------------------------------------------------
[ m_sig,n_k,mu_sig,n_i,n_j ] = Initialize_Parameters();

[ ms,mus,K_y,K_n,ks_i,ks_n ] = Simulate_Data( m_sig,n_k,mu_sig,n_i,n_j );

[ N_iter,i_iter,I_del,M_tol,K_iters ] = Initialize_Estimation_Parameters();

Mu_sig = mu_sig;                            % use actual values (need to estimate these later)
M_sig  = m_sig;

[ Ms,MsAll,MusAll,Mus_ks_n,n_iters ] = Estimate_Parameters(Mu_sig,n_k,M_sig,n_i,...
    K_y,K_n,ks_i,ks_n,N_iter,i_iter,I_del,M_tol,K_iters);

Plot_Results(ms,mus,Ms,MsAll,MusAll,Mus_ks_n,ks_n,n_i,n_j,n_k,i_iter,K_iters,n_iters)


function [ m_sig,n_k,mu_sig,n_i,n_j ] = Initialize_Parameters()
%----------------------------------------------------------------------------------------------
% Initialize parameters
%----------------------------------------------------------------------------------------------
m_sig  = 2;                                 % sigma for the prior of M
n_k    = 500;                               % number of ideas
mu_sig = 1.0;                               % sigma for individual mus
n_i    = 1000/2;                            % number of individuals
n_j    = 50;                                % number of comparisons
n_j_max= n_k*(n_k-1)/2;
n_j    = min(n_j,n_j_max);                  % make sure it does not exceed the maximum possible
fprintf('n_j = %i, n_j (max) = %i\n',n_j,n_j_max)


function [ ms,mus,K_y,K_n,ks_i,ks_n ] = Simulate_Data( m_sig,n_k,mu_sig,n_i,n_j )
%----------------------------------------------------------------------------------------------
% Initialize arrays
%----------------------------------------------------------------------------------------------
ms     = randn(1,n_k)*m_sig;                % mean values
mus    = randn(n_i,n_k)*mu_sig + ms;        % individual mus
               
k_y    = zeros(n_i,n_j);                    % ideas that are preferred
k_n    = k_y;                               %                not preferred
K_y    = k_y;
K_n    = k_n;
ks_i   = cell(n_i,1);
ks_n   = zeros(n_i,1);

%----------------------------------------------------------------------------------------------
% Generate simulated data
%----------------------------------------------------------------------------------------------
for i=1:n_i                                 % for each individual
    j = 0;
    while j<n_j                             % for each comparison
        k_1  = randi(n_k,1);                % first idea
        j_2n = (k_y(i,:)==k_1);             % comparisons where the first idea is selected
        j_2y = (k_n(i,:)==k_1);             %                                     not selected
        k_2n = k_n(i,j_2n);                 % corresponding k_n
        k_2y = k_y(i,j_2y);                 %               k_y
        k_2s = setdiff(1:n_k,[k_2n k_2y k_1]);  % all remaining values
        l_2  = length(k_2s);                    % number of remaining values
        if l_2                                  % if at least one
            k_2  = k_2s(randi(l_2,1));          % second idea
            j          = j + 1;                                      % increment index
            choose_k_1 = rand<normcdf(mus(i,k_1)-mus(i,k_2));        % preference for k=1
            if choose_k_1                   % if k=1 is preferred
                k_y(i,j) = k_1;             % set the indices for this preference
                k_n(i,j) = k_2;
            else                            % otherwise k=2 is preferred
                k_y(i,j) = k_2;
                k_n(i,j) = k_1;
            end
        end
    end
    ks_i{i} = union(k_y(i,:),k_n(i,:));     % get union
    ks_n(i) = length(ks_i{i});              % get number
    [~,K_y(i,:)] = ismember(k_y(i,:),ks_i{i});
    [~,K_n(i,:)] = ismember(k_n(i,:),ks_i{i});
end

flag_print_k = 0;
n_j_max      = 20;

if n_j>n_j_max
    flag_print_k = 0;
    fprintf('\nn_j = %i > %i, Cannot print\n',n_j,n_j_max)
end

if flag_print_k
    fmt = '%2i,%2i';
    for i=1:n_i                                     % for each individual
        fprintf('%3i,',i)
        for j=1:n_j-1                               % for each comparison
            fprintf([fmt ','],k_y(i,j),k_n(i,j))
        end
        fprintf([fmt '\n'],k_y(i,n_j),k_n(i,n_j))
    end
end


function [ N_iter,i_iter,I_del,M_tol,K_iters ] = Initialize_Estimation_Parameters()
%----------------------------------------------------------------------------------------------
% Parameters for estimation algorithm
%----------------------------------------------------------------------------------------------

N_iter = 500/10;                               % number of iterations
i_iter = 100;

I_del  = 0.10;
M_tol  = 0.001;

K_iters= [ 1 30 ];


function [ Ms,MsAll,MusAll,Mus_ks_n,n_iters ] = Estimate_Parameters(Mu_sig,n_k,M_sig,n_i,...
    K_y,K_n,ks_i,ks_n,N_iter,i_iter,I_del,M_tol,K_iters)
%----------------------------------------------------------------------------------------------
% Estimate parameters
%----------------------------------------------------------------------------------------------

Ms     = zeros(1,n_k);                      % initialize estimates for ms
Mus    = cell(n_i,1);
n_iters= length(K_iters);
MusAll = cell(n_iters,i_iter,n_i);
MsAll  = zeros(N_iter,n_k);
Mus_ks_n = nan(n_i,n_k);

for i=1:n_i
    Mus{i} = zeros(1,ks_n(i));
    for iK_iters=1:n_iters
        for I_iter=1:i_iter
            MusAll{iK_iters,I_iter,i} = Mus{i};
        end
    end
end

flag_print = 1;
n_k_max    = 30;
n_i_max    = 100;

if flag_print && (n_k>n_k_max || n_i>n_i_max)
    flag_print = 0;
    fprintf('\nn_k = %i > %i or n_i = %i > %i, cannot print\n',n_k,n_k_max,n_i,n_i_max)
end

if flag_print
    fmt    = '%4i';
    for k=1:n_k, fmt = [fmt ' %7.1f']; end  % format for outputting results to check convergence
    fmt    = [fmt '\n'];
end

%----------------------------------------------------------------------------------------------
% Estimation analysis explained in detail
%----------------------------------------------------------------------------------------------
% K_iter    = counter for the number of global iterations
% Ms(k)     = estimate for the global mean values for idea k for all individuals
% ks_i{i}   = set of ideas that were examined by i
% Msi       = estimate for the global mean values of the set of ideas examined by individual i
% Mu_i(k)   = estimate for the individual mean values for idea k, individual i
% ky(l)     = for comparison l, the idea that was preferred
% kn(l)     = for comparison l, the idea that was not preferred
% I_iter    = counter for the number of iterations for individual i
% i_iter    = maximum number of interations for each individual i
% PSI1      = sum of I*PHI - (I*PHI)' where PSI = phi/PHI (see handwritten notes for explanation)
% Mu_i(new) = Msi + PSI1*Mu_sig^2 (see handwritten notes)
%             this is solved iteratively by using the weighted average:
%             Mu_i = (1-I_del)*Mu_i(old) + I_del*Mu_i(new)
%----------------------------------------------------------------------------------------------

K_iter = 0;

while K_iter<N_iter                         % for each iteration
    K_iter          = K_iter + 1;
    MsAll(K_iter,:) = Ms;                   % current estimates for global means
    Ms_RMS = sqrt(sum(Ms.^2)) + M_tol;      % size of vector 
    for i=1:n_i                             % for each individual
        Mu_i = Mus{i};                      % current estimates for means mu_i
        Msi  = Ms(ks_i{i});                 % global means of ideas considered by i
        ky   = K_y(i,:);                    % preferred ideas (currently same number of comparisons for each i)
        kn   = K_n(i,:);                    % not preferred
        I_iter = 0;
        while I_iter<i_iter                             % for each inner iteration
            I_iter   = I_iter + 1;
            iK_iters = find(K_iter==K_iters);           % check if it is one of K_iters for saving results
            if ~isempty(iK_iters)                       % if it is
                MusAll{iK_iters,I_iter,i}(:) = Mu_i;    %   then save
            end
            PSI1     = calc_PSI1( ky, kn, Mu_i, ks_n(i) ); % sum of I*PSI - (I*PSI)' for each idea considered by i
            Mu_i_new = Msi + PSI1'*Mu_sig^2;            % updated Mu_i
            Mu_del   = Mu_i_new - Mu_i;                 % change in Mu_i  
            Mu_i     = Mu_i + I_del*Mu_del;             % update at reduced increment
            Mu_diff  = sqrt(sum(Mu_del.^2))/Ms_RMS;     % normalized difference
            if Mu_diff < M_tol || isnan(Mu_diff), break, end % if change is small, end
        end        
        Mus{i}(:) = Mu_i;                               % save the result for individual i
        Mus_ks_n(i,ks_i{i}) = Mu_i;                     % save in the full matrix of individual means
    end
    
    Ms_new = nanmean(Mus_ks_n)/(1 + Mu_sig^2/(n_i*M_sig^2)); % average of individual used for the updated estimate for Ms
    Ms_del = Ms_new - Ms;                                    % difference of estimates
    Ms     = Ms_new;                                         % updated estimate
    if sqrt(sum(Ms_del.^2))/Ms_RMS < M_tol, break, end       % if change is small, end
    if flag_print, fprintf(fmt,K_iter,Ms), end
end


function Plot_Results(ms,mus,Ms,MsAll,MusAll,Mus_ks_n,ks_n,n_i,n_j,n_k,i_iter,K_iters,n_iters)
%----------------------------------------------------------------------------------------------
% Plot
%----------------------------------------------------------------------------------------------
n_param = sum(~isnan(Mus_ks_n(:)));

fig = figure(102); fig.Name = 'Comparison'; fig.Position = [500 180 1000 800]; fig.Color = 'w';
subplot(1,2,1), plot(ms, Ms, 'o')
xlabel(texlabel('m_k (actual)')), ylabel(texlabel('m_k (MAP)'))
title(sprintf('MAP results\nn_i = %i individuals, n_k = %i ideas, n_j = %i comparisons\nPrior means',...
    n_i,n_k,n_j))
subplot(1,2,2), plot(mus,Mus_ks_n,'.')
xlabel(texlabel('mu_k^i (actual)')), ylabel(texlabel('mu_k^i (MAP)'))
title(sprintf('Prior individual means\nnumber of parameters = %i',n_param))

np  = n_iters+1;                                   % number of plots
fig = figure(110); fig.Name = 'Convergence';
subplot(np,1,1), plot(MsAll), xlabel('iteration'), ylabel('m_k'), title('Convergence of m_k')
i   = 2;                                           % individual of interest
for iK_iters=1:n_iters
    K_iter = K_iters(iK_iters);
    subplot(np,1,1+iK_iters)
    plot(reshape(cell2mat(MusAll(iK_iters,:,i)),ks_n(i),i_iter)')
    xlabel('iteration'), ylabel(texlabel('mu_k^i'))
    title(texlabel(sprintf('Convergence of mu_k^i, at K_iter = %i',K_iter)))
end

fprintf('\nNumber of estimated parameters: m_k (%i), m_k^i (%i)\n',...
    n_k,n_param)


function PSI1 = calc_PSI1( k_y, k_n, Mu, n_k )
%----------------------------------------------------------------------------------------------
% Calculate psi1
%----------------------------------------------------------------------------------------------
PSI   = zeros(n_k);                           % initialize
kyn   = sub2ind([n_k n_k],k_y,k_n);           % get index
mu_D  = Mu(k_y) - Mu(k_n);                    % difference in mus
PSI(kyn) = normpdf(mu_D)./normcdf(mu_D);      % set the matrix value

PSI1  = sum(PSI-PSI',2);                      % sum of I*PSI - (I*PSI)'
