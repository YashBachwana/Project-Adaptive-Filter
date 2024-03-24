clc ; clear all ; close all ;
input_length = 10000 ;
SNR = 30 ;
for iter = 1 : 1000
    iter
input = rand(1,input_length) - 0.5;
input = 2 * input;

% Calculate the mean of the input
mean_input = mean(input);

% Adjust the mean to 0
input = input - mean_input;
system_noise = awgn(input,SNR)-input ;

% system_output = g(input) ;
numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
denominator_coefficients = [1, -1.99, 1.572, -0.4583];
system_output = sin(filter(numerator_coefficients, denominator_coefficients, input)) + system_noise ;

x_ = -2.5 : 0.05 : 2.5 ;
x_out = sin(filter(numerator_coefficients, denominator_coefficients, x_)) ;

x1 = -2.5:0.05:-1.5 ; x2 = -2:0.05:-1 ; x3 = -1.5:0.05:-0.5 ; x4 = -1:0.05:0 ; x5 = -0.5:0.05:0.5 ;
x6 = 0:0.05:1 ; x7 = 0.5:0.05:1.5 ; x8 = 1:0.05:2 ; x9 = 1.5:0.05:2.5 ;
x = {x1,x2,x3,x4,x5,x6,x7,x8,x9} ;

Y = x ; Z = x ;

Q = (101:-20:-72)';
P = 2 ; del_x = 0.05 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;
mu = 0.025 ; n = 9 ; 

for i = 1 : input_length
    u = (input(i)/del_x) - floor(input(i)/del_x);
    parameter_array_u = [(u)^2,u,1];
    
    for j = 1 : n
        index = floor(input(i)/del_x) + (Q(j) - 1)/2 + 1 ;
        if (index >= 1 && index <= 18)
            out = parameter_array_u * C * Y{j}(index : index + 2)' ;
            error = system_output(i) - out ;
            if (index >= 6)
                err(i) = error ;
            end 
            Y{j}(index:index + 2) = Y{j}(index : index + 2)' + mu * error * C' * parameter_array_u' ;
            
        end 
    end 
end 
err_ensemble_kernel(iter,:) = err .^ 2 ;
end 
plot(10 * log10(mean(err_ensemble_kernel)),'b') ;

x = -2.5:0.05:2.5 ;
z = sin(filter(numerator_coefficients, denominator_coefficients, x)) ;
y = zeros(1,101) ; j = 6 ;

k = 1 ; j = 1 ;
var = 0 ;
for i = 1 : 101
    var = i ;
    if (k < 15)
        y(i) = Y{j}(k) ;
        k = k + 1 ;
    elseif k == 15
        y(i) = Y{j}(k) ;
        j = j + 1 ;
        if (j > 9)
            break 
        end 
        k = 6 ;    
    end  
end 
y(var + 1:101) = Y{n}(16:21) ;

