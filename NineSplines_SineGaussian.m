clc ; clear all ; close all ;
input_length = 10000 ;
SNR = 30 ;

for iter = 1 : 1000
    iter
input = rand(1,input_length) - 0.5;
input = 2 * input;

mean_input = mean(input);

input = input - mean_input;
system_noise = awgn(input,SNR)-input ;

system_output = g(input) + system_noise;
x_ = -2.5 : 0.05 : 2.5 ;

x1 = -2.5:0.05:-1.5 ; x2 = -2:0.05:-1 ; x3 = -1.5:0.05:-0.5 ; x4 = -1:0.05:0 ; x5 = -0.5:0.05:0.5 ;
x6 = 0:0.05:1 ; x7 = 0.5:0.05:1.5 ; x8 = 1:0.05:2 ; x9 = 1.5:0.05:2.5 ;
x = {x1,x2,x3,x4,x5,x6,x7,x8,x9} ;

Y = x ; Z = x ;

Q = (101:-20:-72)';
P = 2 ; del_x = 0.05 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;
mu = 0.25 ; n = 9 ; 

    for i = 1 : input_length
        i   ; 
        u = (input(i)/del_x) - floor(input(i)/del_x);
        parameter_array_u = [(u)^2,u,1];
        
        for j = 1 : n
            index = floor(input(i)/del_x) + (Q(j) - 1)/2 + 1 ;
            if (index >= 1 && index <= 18)
                out = parameter_array_u * C * Y{j}(index : index + 2)' ;
                error = system_output(i) - out ;
                
                if (index >= 6)
                    err(i) = error ;
                    Y{j}(index:index + 2) = Y{j}(index : index + 2)' + mu * error * C' * parameter_array_u' ; 
                    break ;
                end 
                
                
            end 
        end 
    
    end 
    err_ensemble_kernel(iter,:) = err .^ 2 ;
end 
plot(10 * log10(mean(err_ensemble_kernel)),'b') ;  



function system_output = g(x)
    input_length = length(x) ;
    system_output = zeros(1,input_length) ;
    for i = 1:input_length
        system_output(i) = f(x(i)) ;
    end
    system_output = sin(system_output) ;
end 

function physical_output = f(x)
    k1 = 40 ; k2 = 30 ; h1 = 5 ; h2 = 0.5 ; c1 = -0.8 ; c2 = 0.5 ;
    physical_output = sin(k1 * exp(-1 * ((x - c1) ^ 2) / (2 * h2^2)) +  k2 * exp(-1 * ((x - c2) ^ 2) / (2 * h1^2)));
end 