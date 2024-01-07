clc ; clear all ; close all ;

for iter = 1 : 500
    iter

rand_val = rand(1,1000) - 0.5 ;

neeta = 1.25 ; 
phy_arr_in = [0,0,0,0] ;
phy_array_out = [0,0,0,0] ;
d(1) = 0 ; h = 10000000 ; e(1) = [0] ; e_linear(1) = [0] ;

%Linear Filter Initialisation 
W = [0,0,0] ;
X_phy = [0,0,0] ;
phy_arr = [1,2,1] ;

% kernel initialisations
a(1) = neeta * d(1) ;
C(1) = {phy_arr_in} ;

noise = awgn(rand_val,30) - rand_val ;

for i = 2 : length(rand_val)
  
    phy_arr_in(i + 3) = rand_val(i) ;
    if i <= 0 
        phy_out = X_phy * phy_arr' ;
        phy_array_out(i + 3) = phy_out ; 
        phy_out = phy_out + noise(i);
    
    else 
    c1 = 0.0154 ; c2 = 0.0462 ; c3 = c1 ; c4 = 1.99 ; c5 = 1.572 ; c6 = 0.4583 ;  
    phy_out = c1 * phy_arr_in(i + 3) + c2 * phy_arr_in(i + 2) + c2 * phy_arr_in(i + 1) + c3 * phy_arr_in(i) + c4 * (phy_array_out(i + 2)) - c5 * phy_array_out(i + 1) + c6 * phy_array_out(i) ;
    phy_array_out(i + 3) = phy_out ; 
    phy_out = sin(phy_out);
    phy_out = phy_out + noise(i);
    end 
    sys_out = 0 ;
    for j = 1 : i - 1
        sys_out = sys_out + a(j) * exp(-1 * ((rand_val(i) - rand_val(j)) ^ 2) / (2 * h^2) )  ;
    end 
    

    X_phy = [rand_val(i), X_phy(1:length(X_phy) - 1)] ;
    sys_out_linear = X_phy * W' ;
    e_linear(i) = phy_out - sys_out_linear ;
    W = W + 2 * 0.005 * e_linear(i) * X_phy ;
    
   
    e(i) = phy_out - sys_out ;
   
    a(i) = neeta * e(i) ; 
end 
err_ensemble_kernel(iter,:) = e .^ 2 ;
err_ensemble_linear(iter,:) = e_linear .^ 2 ;
end 

plot(10 * log10(mean(err_ensemble_kernel)),'r') ;
hold on ;
plot(10 * log10(mean(err_ensemble_linear)),'k') ;
hold off ; 
legend('kernel', 'linear', 'Location', 'southoutside', 'Orientation', 'horizontal');
