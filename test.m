
for i = 1:size(V,2)-1
diffV(:,i) = vV(:,i)-vV(:,i+1);
end
for i =1:size(diffV, 2)
diffnormv(i) = norm(diffV(:,i));
end

for i = 1:size(V,2)-1
diffV(:,i) = vV_mp(:,i)-vV_mp(:,i+1);
end
for i =1:size(diffV, 2)
diffnormvm(i) = norm(diffV(:,i));
end

plot(diffnormvm-diffnormv)


%%

mpc = case9;
[V, V_mp, distance1] = socp_experiment(mpc, 0);
[V, V_mp, distance] = socp_experiment(mpc, 0, V_mp);
% Vhat=V_mp;

%%


sci = 1; % intensity of wrong measurements
s = 0.00; %0.01; % Noise level
op = 0; % 0.1;
% M = getM(xhat)


format long
warning off

addpath('.\utils');

alpha = 0; 

nb = size(mpc.bus, 1);

busL = mpc.bus(:,1);
busN = sparse(busL, ones(nb,1), (1:nb));

mpc.bus(:,1) = busN(mpc.bus(:,1));
mpc.branch = mpc.branch(find(mpc.branch(:,11)),:);
mpc.branch(:,1) = busN(mpc.branch(:,1));
mpc.branch(:,2) = busN(mpc.branch(:,2));
mpc.gen(:,1) = busN(mpc.gen(:,1));

nl = size(mpc.branch, 1);
ng = size(mpc.gen, 1);

[Y,Yf,Yt] = makeYbus(mpc);

Cf = sparse(1:nl, mpc.branch(:,1), 1, nl, nb);
Ct = sparse(1:nl, mpc.branch(:,2), 1, nl, nb);
Cg = sparse(1:ng, mpc.gen(:,1), 1, ng, nb);

ff = mpc.branch(:,1);
tt = mpc.branch(:,2);

[tw, perm, bags] = Permutation(1:nb,ff,tt,alpha);
nc = size(bags,1);
%%%

bags = cell(1);
for kk = 1 : nl
    bags{kk} = [ff(kk),tt(kk)];
end
nc = nl;

ndx = sparse(nb,nb);
for kk = 1 : nc
    ndx(bags{kk},bags{kk}) = 1;
end

nn = sum(sum(ndx));

B = sparse(nn,nc);
bb = zeros(nc,1);
for kk = 1 : nc
%     [kk,1]
    
    bb(kk) = length(bags{kk});
    cc = sparse(nb,nb);
    cc(bags{kk},bags{kk}) = v2m(ones(bb(kk)^2,1),ones(bb(kk)));
    B(:,kk) = m2v(cc,ndx);
end

%%%%%% result of matpower
res_mp = runpf(mpc);
V_mp   = res_mp.bus(:,8) .* exp((res_mp.bus(:,9) - res_mp.bus(find(mpc.bus(:,2) == 3),9))*pi*1i/180);
pf_mp = res_mp.branch(:,14)/mpc.baseMVA;
qf_mp = res_mp.branch(:,15)/mpc.baseMVA;
pt_mp = res_mp.branch(:,16)/mpc.baseMVA;
qt_mp = res_mp.branch(:,17)/mpc.baseMVA;
W_mp  = V_mp*V_mp';
x_mp  = m2v(W_mp,ndx);
%%%%%%

Yp = cell(nb,1);
Yq = cell(nb,1);
Yfp = cell(nl,1);
Yfq = cell(nl,1);
Ytp = cell(nl,1);
Ytq = cell(nl,1);
for ii = 1 : nb
%     [ii,2]
    
    en     = sparse(nb,1);
    en(ii) = 1;
    
    Yp{ii} = (Y'*en*en' + en*en'*Y)/2;
    Yq{ii} = (Y'*en*en' - en*en'*Y)/2i;
end

%%%%%
for ll = 1 : nl
%     [ll,3]
    
    ef = sparse(nb,1);
    ef(ff(ll)) = 1;
    
    et = sparse(nb,1);
    et(tt(ll)) = 1;

    el     = sparse(nl,1);
    el(ll) = 1;
    
    Yfp{ll} = (Yf'*el*ef' + ef*el'*Yf)/2;
    Yfq{ll} = (Yf'*el*ef' - ef*el'*Yf)/2i;
    
    Ytp{ll} = (Yt'*el*et' + et*el'*Yt)/2;
    Ytq{ll} = (Yt'*el*et' - et*el'*Yt)/2i;
end

Dp = sparse(nb, nn);
Dq = sparse(nb, nn);
Dv = sparse(nb, nn);
for ii = 1 : nb
%     [ii,4]
    
    ee = sparse(nb,nb);
    ee(ii,ii) = 1;
    Dv(ii,:) = m2v(ee,ndx)';
    Dp(ii,:) = m2v(Yp{ii},ndx)';
    Dq(ii,:) = m2v(Yq{ii},ndx)';
end

Dfp = sparse(nl, nn);
Dfq = sparse(nl, nn);
Dtp = sparse(nl, nn);
Dtq = sparse(nl, nn);
for ll = 1 : nl
%     [ll,5]
    
    Dfp(ll,:) = m2v(Yfp{ll},ndx)';
    Dfq(ll,:) = m2v(Yfq{ll},ndx)';
    Dtp(ll,:) = m2v(Ytp{ll},ndx)';
    Dtq(ll,:) = m2v(Ytq{ll},ndx)';
end

dv  = Dv * x_mp;
dp  = Dp * x_mp;
dq  = Dq * x_mp;
dfp = Dfp * x_mp;
dfq = Dfq * x_mp;
dtp = Dtp * x_mp;
dtq = Dtq * x_mp;


ff_x = zeros(nl,1);
tt_x = zeros(nl,1);
ft_xr = zeros(nl,1);
ft_xi = zeros(nl,1);
for ll = 1 : nl
%     [ll,6]
    
    ef = sparse(nb,nb);
    et = sparse(nb,nb);
    fti = sparse(nb,nb);
    ftr = sparse(nb,nb);
    
    ef(ff(ll),ff(ll)) = 1;
    et(tt(ll),tt(ll)) = 1;
    ftr(ff(ll),tt(ll)) = 1;
    ftr(tt(ll),ff(ll)) = 1;
    fti(ff(ll),tt(ll)) = +1i;
    fti(tt(ll),ff(ll)) = -1i;
    
    ff_x(ll) = find(m2v(ef,ndx));
    tt_x(ll) = find(m2v(et,ndx));
    ft_xr(ll) = find(m2v(ftr,ndx));
    ft_xi(ll) = find(m2v(fti,ndx));
end

%%%%%%%%%%

vv = [1:nb];

pp  = [1:nb];

qq  = [1:nb];

ppf = [1:nl];
qqf = [1:nl];
ppt = [1:nl];
qqt = [1:nl];

oo = length(vv) + length(pp) + length(qq) + ...
     length(ppf) + length(qqf) + length(ppt) + length(qqt);

DD = [Dv(vv,:); Dp(pp,:); Dq(qq,:); Dfp(ppf,:); Dfq(qqf,:); Dtp(ppt,:); Dtq(qqt,:)];
dd = [dv(vv); dp(pp); dq(qq); dfp(ppf); dfq(qqf); dtp(ppt); dtq(qqt);];

%%%%%%%%%% 

% generate dense noise
    
ssv  = s*abs(Dv*x_mp);
ssp  = s*abs(Dp*x_mp);
ssq  = s*abs(Dq*x_mp);
ssfp = s*abs(Dfp*x_mp);
ssfq = s*abs(Dfq*x_mp);
sstp = s*abs(Dtp*x_mp);
sstq = s*abs(Dtq*x_mp);

sv  = ssv.*randn(nb,1);
sp  = ssp.*randn(nb,1);
sq  = ssq.*randn(nb,1);
sfp = ssfp.*randn(nl,1);
sfq = ssfq.*randn(nl,1);
stp = sstp.*randn(nl,1);
stq = sstq.*randn(nl,1);

ss = [sv(vv); sp(pp); sq(qq); sfp(ppf); sfq(qqf); stp(ppt); stq(qqt);];

sc = zeros(oo,1);
sc(randsample(oo,floor(op*oo))) = 1; % portion of wrong measurements

%%%%%%%%%%

ssb = ss + sci*sc; % noise = dense noise + sparse noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (nargin == 2)
%         aa = max(eigs(imag(Y),nb));
%         M  = -imag(Y)+sparse(1:nb,1:nb,aa,nb,nb);
% elseif (nargin == 3)
%         M = getM(Vhat, ndx);
% end

M = getM(Vhat, ndx);
% M = -triu(ndx,1);
% M = (M + M')/2;

% M = eye(nb);

% Obtain W

Dm = m2v(M,ndx)';
M
Dm
v2m(Dm,ndx)-M

%%

n = length(ndx);
    l = length(V);
    d = trace(ndx);
    M = sparse(n,n);
    
    [fxd,txd] = find(triu(ndx));
    [fx,tx] = find(triu(ndx,1));
    
%     M = M + sparse(fxd, txd, V(1:(l+n)/2), n, n)/sqrt(2);
%     M = M + sparse(txd, fxd, V(1:(l+n)/2), n, n)/sqrt(2);
%     M = M + 1i*sparse(fx, tx, V(1+(l+n)/2:l), n, n)/sqrt(2);
%     M = M - 1i*sparse(tx, fx, V(1+(l+n)/2:l), n, n)/sqrt(2);
    M = M + sparse(fxd, txd, V(1:(l+d)/2), n, n)/sqrt(2);
    M = M + sparse(txd, fxd, V(1:(l+d)/2), n, n)/sqrt(2);
    M = M + 1i*sparse(fx, tx, V(1+(l+d)/2:l), n, n)/sqrt(2);
    M = M - 1i*sparse(tx, fx, V(1+(l+d)/2:l), n, n)/sqrt(2);
    M = M + M .* sparse(1:n,1:n,1/sqrt(2)-1,n,n);

%%
    
n = length(ndx);
MM = M - M.*sparse(1:n,1:n,1-1/sqrt(2),n,n);
V = sqrt(2) * [real(MM(triu(ndx)~=0)); imag(MM(triu(ndx,1)~=0))];

%%

% % mpc = case1354pegase;
% mpc = case39;
% rng(42)
% 
% format long
% warning off
% 
% addpath('.\utils');
% 
% alpha = 0; 
% 
% nb = size(mpc.bus, 1);
% 
% busL = mpc.bus(:,1);
% busN = sparse(busL, ones(nb,1), (1:nb));
% 
% mpc.bus(:,1) = busN(mpc.bus(:,1));
% mpc.branch = mpc.branch(find(mpc.branch(:,11)),:);
% mpc.branch(:,1) = busN(mpc.branch(:,1));
% mpc.branch(:,2) = busN(mpc.branch(:,2));
% mpc.gen(:,1) = busN(mpc.gen(:,1));
% 
% nl = size(mpc.branch, 1);
% ng = size(mpc.gen, 1);
% 
% [Y,Yf,Yt] = makeYbus(mpc);
% 
% Cf = sparse(1:nl, mpc.branch(:,1), 1, nl, nb);
% Ct = sparse(1:nl, mpc.branch(:,2), 1, nl, nb);
% Cg = sparse(1:ng, mpc.gen(:,1), 1, ng, nb);
% 
% ff = mpc.branch(:,1);
% tt = mpc.branch(:,2);
% 
% [tw, perm, bags] = Permutation(1:nb,ff,tt,alpha);
% nc = size(bags,1);
% %%%
% 
% bags = cell(1);
% for kk = 1 : nl
%     bags{kk} = [ff(kk),tt(kk)];
% end
% nc = nl;
% 
% ndx = sparse(nb,nb);
% for kk = 1 : nc
%     ndx(bags{kk},bags{kk}) = 1;
% end
% 
% % ndx = ones(nb,nb);
% 
% nn = sum(sum(ndx));
% 
% B = sparse(nn,nc);
% bb = zeros(nc,1);
% for kk = 1 : nc
% %     [kk,1]
%     
%     bb(kk) = length(bags{kk});
%     cc = sparse(nb,nb);
%     cc(bags{kk},bags{kk}) = v2m(ones(bb(kk)^2,1),ones(bb(kk)));
%     B(:,kk) = m2v(cc,ndx);
% end
% 
% %%%%%% result of matpower
% res_mp = runpf(mpc);
% V_mp   = res_mp.bus(:,8) .* exp((res_mp.bus(:,9) - res_mp.bus(find(mpc.bus(:,2) == 3),9))*pi*1i/180);
% pf_mp = res_mp.branch(:,14)/mpc.baseMVA;
% qf_mp = res_mp.branch(:,15)/mpc.baseMVA;
% pt_mp = res_mp.branch(:,16)/mpc.baseMVA;
% qt_mp = res_mp.branch(:,17)/mpc.baseMVA;
% W_mp  = V_mp*V_mp';
% x_mp  = m2v(W_mp,ndx);
% %%%%%%
% 
% Yp = cell(nb,1);
% Yq = cell(nb,1);
% Yfp = cell(nl,1);
% Yfq = cell(nl,1);
% Ytp = cell(nl,1);
% Ytq = cell(nl,1);
% for ii = 1 : nb
% %     [ii,2]
%     
%     en     = sparse(nb,1);
%     en(ii) = 1;
%     
%     Yp{ii} = (Y'*en*en' + en*en'*Y)/2;
%     Yq{ii} = (Y'*en*en' - en*en'*Y)/2i;
% end
% 
% %%%%%
% for ll = 1 : nl
% %     [ll,3]
%     
%     ef = sparse(nb,1);
%     ef(ff(ll)) = 1;
%     
%     et = sparse(nb,1);
%     et(tt(ll)) = 1;
% 
%     el     = sparse(nl,1);
%     el(ll) = 1;
%     
%     Yfp{ll} = (Yf'*el*ef' + ef*el'*Yf)/2;
%     Yfq{ll} = (Yf'*el*ef' - ef*el'*Yf)/2i;
%     
%     Ytp{ll} = (Yt'*el*et' + et*el'*Yt)/2;
%     Ytq{ll} = (Yt'*el*et' - et*el'*Yt)/2i;
% end
% 
% Dp = sparse(nb, nn);
% Dq = sparse(nb, nn);
% Dv = sparse(nb, nn);
% for ii = 1 : nb
% %     [ii,4]
%     
%     ee = sparse(nb,nb);
%     ee(ii,ii) = 1;
%     Dv(ii,:) = m2v(ee,ndx)';
%     Dp(ii,:) = m2v(Yp{ii},ndx)';
%     Dq(ii,:) = m2v(Yq{ii},ndx)';
% end
% 
% Dfp = sparse(nl, nn);
% Dfq = sparse(nl, nn);
% Dtp = sparse(nl, nn);
% Dtq = sparse(nl, nn);
% for ll = 1 : nl
% %     [ll,5]
%     
%     Dfp(ll,:) = m2v(Yfp{ll},ndx)';
%     Dfq(ll,:) = m2v(Yfq{ll},ndx)';
%     Dtp(ll,:) = m2v(Ytp{ll},ndx)';
%     Dtq(ll,:) = m2v(Ytq{ll},ndx)';
% end
% 
% dv  = Dv * x_mp;
% dp  = Dp * x_mp;
% dq  = Dq * x_mp;
% dfp = Dfp * x_mp;
% dfq = Dfq * x_mp;
% dtp = Dtp * x_mp;
% dtq = Dtq * x_mp;
% 
% 
% ff_x = zeros(nl,1);
% tt_x = zeros(nl,1);
% ft_xr = zeros(nl,1);
% ft_xi = zeros(nl,1);
% for ll = 1 : nl
% %     [ll,6]
%     
%     ef = sparse(nb,nb);
%     et = sparse(nb,nb);
%     fti = sparse(nb,nb);
%     ftr = sparse(nb,nb);
%     
%     ef(ff(ll),ff(ll)) = 1;
%     et(tt(ll),tt(ll)) = 1;
%     ftr(ff(ll),tt(ll)) = 1;
%     ftr(tt(ll),ff(ll)) = 1;
%     fti(ff(ll),tt(ll)) = +1i;
%     fti(tt(ll),ff(ll)) = -1i;
%     
%     ff_x(ll) = find(m2v(ef,ndx));
%     tt_x(ll) = find(m2v(et,ndx));
%     ft_xr(ll) = find(m2v(ftr,ndx));
%     ft_xi(ll) = find(m2v(fti,ndx));
% end
% 
% %%%%%%%%%%
% 
% vv = [1:nb];
% 
% pp  = [1:nb];
% 
% qq  = [1:nb];
% 
% ppf = [1:nl];
% qqf = [1:nl];
% ppt = [1:nl];
% qqt = [1:nl];
% 
% oo = length(vv) + length(pp) + length(qq) + ...
%      length(ppf) + length(qqf) + length(ppt) + length(qqt);
% 
% DD = [Dv(vv,:); Dp(pp,:); Dq(qq,:); Dfp(ppf,:); Dfq(qqf,:); Dtp(ppt,:); Dtq(qqt,:)];
% dd = [dv(vv); dp(pp); dq(qq); dfp(ppf); dfq(qqf); dtp(ppt); dtq(qqt);];
% 
% %%%%%%%%%% 
% 
% % generate dense noise
% 
% s = 0.00; %0.01; % Noise level
%     
% ssv  = s*abs(Dv*x_mp);
% ssp  = s*abs(Dp*x_mp);
% ssq  = s*abs(Dq*x_mp);
% ssfp = s*abs(Dfp*x_mp);
% ssfq = s*abs(Dfq*x_mp);
% sstp = s*abs(Dtp*x_mp);
% sstq = s*abs(Dtq*x_mp);
% 
% sv  = ssv.*randn(nb,1);
% sp  = ssp.*randn(nb,1);
% sq  = ssq.*randn(nb,1);
% sfp = ssfp.*randn(nl,1);
% sfq = ssfq.*randn(nl,1);
% stp = sstp.*randn(nl,1);
% stq = sstq.*randn(nl,1);
% 
% ss = [sv(vv); sp(pp); sq(qq); sfp(ppf); sfq(qqf); stp(ppt); stq(qqt);];
% 
% op = 0.1;
% sci = 1; % intensity of wrong measurements
% 
% sc = zeros(oo,1);
% sc(randsample(oo,floor(op*oo))) = 1; % portion of wrong measurements
% 
% 
% %%%%%%%%%%
% 
% ssb = ss + sci*sc; % noise = dense noise + sparse noise
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % aa = max(eigs(imag(Y),nb));
% % M  = -imag(Y)+sparse(1:nb,1:nb,aa,nb,nb);
% 
% M = getM(V_mp, ndx);
% 
% 
% Dm = m2v(M,ndx)';
% 
% mu  = 5e-1; %regularization parameter
% cvx_begin
% % cvx_solver_settings('maxit',150,'gaptol',1e-15,'inftol',1e-15,'steptol',1e-15);
% cvx_precision best
% 
% variable x(nn)
% variable h(oo)
% 
% r = dd + ssb - DD * x;
% 
% minimize(Dm*x + mu*sum(abs(h)));
% 
% subject to
%     
%     %WLS
% %     (r.^2) <= h;
% 
%     %WLAV
%     abs(r) <= h;
% 
% %     for kk = 1 : nc-1
% %         v2m(x(find(B(:,kk))),ones(bb(kk))) == hermitian_semidefinite(bb(kk));
% %     end
% 
%     for kk = 1 : nc
%         v2m(x(find(B(:,kk))),ones(bb(kk))) == hermitian_semidefinite(bb(kk));
%     end
%     
% cvx_end
% 
% W = v2m(x,ndx);
% 
% max(abs(x - x_mp))
% 
% cvx_begin
% cvx_precision best
% 
% variable tb(nb)
% expression tl(nl)
% 
%     tl = tb(ff) - tb(tt);
%     cost_lp = sum(sum(abs( tl - angle(diag(W(ff, tt))) )));
% 
% minimize( cost_lp );
% 
% subject to
% 
%     tb(find(mpc.bus(:,2) == 3)) == 0;
% 
% cvx_end
% 
% V = full(sqrt(diag(W)).*exp(1i*tb));
% 
% result = sqrt(sum((abs(V_mp-V)).^2)/nb)
