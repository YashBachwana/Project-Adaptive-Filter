clc;
clear all;
close all;

for iter = 1 : 100
    iter
input_length = 10000 ;
SNR = 30 ;
input = rand(1,input_length) - 0.5;
system_noise = awgn(input,SNR)-input ;
system_output = zeros(1,input_length) ;
numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
denominator_coefficients = [1, -1.99, 1.572, -0.4583];
system_output = sin(filter(numerator_coefficients, denominator_coefficients, input));
n = 100 ; 
Q = 21 ; P = 2 ; del_x = 0.2 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;

u = zeros(1,input_length);

index_array = ones(n, input_length) ; 

multiplier = ones(1,n) ; 

error = zeros(1,input_length);

model_output = zeros(1,input_length);
splines_model_output = zeros(n,input_length) ; 

parameter_array_u = zeros(1,P+1);
parameter_array_u_der = zeros(1,P+1);

control_point_array = cell(n, 1);

a = - 2.5 ;
for i = 1 : n
    control_point_array{i} = (a : 5 / (n * 20) : a + (5 / n))' ; 
    a = a + (5 / n) ; 
end 

for i = 1:input_length
    u(i) = (input(i)/del_x) - floor(input(i)/del_x);

    for j = 1:10
        index_array(j,i) = floor(input(i)/del_x)+(Q-1)/2;
    end 

    parameter_array_u = [(u(i))^2,u(i),1];
    parameter_array_u_der = [2*u(i),1,0];
    
    for j = 1:10
        splines_model_output(j,i) = multiplier(j) * parameter_array_u * C * control_point_array{j}(index_array(j,i):index_array(j,i)+2,1) ;
    end 

    model_output(i) = system_noise(i) ; 
    for j = 1 : 10
        model_output(i) = model_output(i) + splines_model_output(j,i) ; 
    end 

    error(i) = system_output(i) - model_output(i) ;
    
    for j = 1:10
        control_point_array{j}(index_array(j,i):index_array(j,i) + 2,1) = control_point_array{j}(index_array(j,i):index_array(j,i) + 2,1) + 0.025 * multiplier(j) * error(i) * C' * parameter_array_u' ; 
    
    end 
    
    for j = 1:10
        multiplier(j) = multiplier(j) + 0.025 * error(i) * parameter_array_u * C * control_point_array{j}(index_array(j,i):index_array(j,i) + 2,1) ; 
    end 
end
err_ensemble_kernel(iter,:) = error .^ 2 ;
system_out_ensemble(iter,:) = system_output .^ 2 ;
model_out_ensemble(iter,:) = model_output .^ 2 ;
end
plot(10 * log10(mean(err_ensemble_kernel)),'b') ; hold on ; 



for iter = 1 : 100
    iter
input_length =10000 ;
SNR = 30 ;
input = rand(1,input_length) - 0.5;
system_noise = awgn(input,SNR)-input ;
system_output = zeros(1,input_length) ;
numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
denominator_coefficients = [1, -1.99, 1.572, -0.4583];
system_output = filter(numerator_coefficients, denominator_coefficients, input);

error = zeros(1,input_length);

model_w = [1,zeros(1,3)];
model_tap = zeros(1,length(model_w));
model_output = zeros(1,input_length);

for i = 1:input_length

    model_tap = [input(i) model_tap(1:end-1)];
    model_output(i) = model_tap * model_w';

    error(i) = system_output(i) - model_output(i) ;

    model_w = model_w + 0.025 * error(i) * model_tap ; 
end
err_ensemble_kernel(iter,:) = error .^ 2 ;
system_out_ensemble(iter,:) = system_output .^ 2 ;
model_out_ensemble(iter,:) = model_output .^ 2 ;
end

plot(10 * log10(mean(err_ensemble_kernel)),'r') ; hold on ; 


legend('Spline System','Linear System')