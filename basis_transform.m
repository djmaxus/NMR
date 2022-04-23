%% Alex's suggesstion on source-term displacement
%% Proposition
%% 
% * Principal components $\Phi$ can be decomposed over fundamental solution 
% $G$, which corresponds to numerical solution $U=G$: $\Phi = G D$
% * Slightly displaced source term corresponds to $\hat{G} = G + \Delta G$
% * *The suggestion is* to calculate displaced coefficients $\hat{D} = D + \Delta 
% D$ and then displaced principal components $\hat{\Phi} = \hat{G} \hat{D}$
% * The declared advantage is a presence of physically-accurate nonlinearities, 
% camptured from $U$, within $D$ and, so, within $\hat{D}$
%% Derivations
% Common
% Principal components and SVD
% $$U = \Phi \Sigma \Upsilon^T$$
% 
% $$K_\text{f} = U U^T$$
% 
% $$K_\text{f} \Phi =  \Phi \Sigma^2$$
% Assumption
% $$\Phi = U D$$
% How we approached it previously
% $$K_\text{f} \times U D = U (U^T U D)= U (D \Sigma^2)$$
% 
% $$U^T U D = D \Sigma^2$$
% First-order perturbation
% $$\hat{U} = U + \Delta U$$
% 
% $$\hat{D} = D + \Delta D$$
% 
% $$V := U^T U$$
% 
% $$\hat{V} = V + \Delta V$$
% First-order solution
% $$(V + \Delta V) (D + \Delta D) = (D + \Delta D) \Sigma^2$$
% 
% $$V \Delta D + D \Delta V = \Delta D \Sigma^2 \\(V - \Sigma^2 ) \mathbf{\Delta 
% D} = D \Delta V$$
% How I see it after one careful notice
% $$K_\text{f} \times U D = U (U^T U D)= U (D \Sigma^2)$$
% 
% $U^T U D = D \Sigma^2$ $\Rightarrow$$D = \Upsilon$ (right singular vectors)
% 
% $$\hat{D} = \hat{\Upsilon} : \hat{V} = \hat{\Upsilon} \hat{\Sigma}^2 \hat{\Upsilon}^T$$
%% Analysis
% "First-order solution" just neglects one term of SVD equality and also a change 
% of singular values $\Sigma^2$. This is unneccesary simplification since SVD 
% relation is some even more accurate equation with known solution. The declared 
% preservation of physically-correct nonlinearities seems making no sense from 
% now on, because it only depends on physical correctness of snapshots in $\hat{V}$.
% 
% Consequently, the more crucial issue is still to find such$\Delta G$ that 
% would be accurate enough to represent a system with displaced source terms. 
% That is exactly what I was focusing on and prefer to carry on.
% 
%