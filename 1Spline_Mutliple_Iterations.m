clc;
clear all;
close all;

for iter = 1 : 50
    iter
input_length = 50000 ;
SNR = 30 ;
input = rand(1,input_length) - 0.5;
system_noise = awgn(input,SNR)-input ;
system_output = zeros(1,input_length) ;
for i = 1:input_length
    system_output(i) = f(input(i)) ;
end

Q = 21 ; P = 2 ; del_x = 0.2 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;

u = zeros(1,input_length);

index_i1 = ones(1,input_length);

multiplier1 = 1 ;

error = zeros(1,input_length);

model_output = zeros(1,input_length);
model_output1 = zeros(1,input_length);

parameter_array_u = zeros(1,P+1);
parameter_array_u_der = zeros(1,P+1);

control_point_array1 = (-2:0.2:2)';

for i = 1:input_length
    u(i) = (input(i)/del_x) - floor(input(i)/del_x);

    index_i1(i) = floor(input(i)/del_x)+(Q-1)/2; 

    parameter_array_u = [(u(i))^2,u(i),1];
    parameter_array_u_der = [2*u(i),1,0];

    model_output1(i) = multiplier1 * parameter_array_u * C * control_point_array1(index_i1(i):index_i1(i)+2,1) ;
    model_output(i) = model_output1(i) + system_noise(i); 

    error(i) = system_output(i) - model_output(i) ;

    control_point_array1(index_i1(i):index_i1(i)+2 , 1) = control_point_array1(index_i1(i):index_i1(i) + 2 , 1) + 0.025 * 2 * multiplier1 * error(i) * C' * parameter_array_u' ;
    multiplier1 = multiplier1 + 0.025 * 2 * error(i) * parameter_array_u * C * control_point_array1(index_i1(i):index_i1(i)+2,1) ;
end
% 
% plot(system_output,'k');
% hold on;
% plot(model_output,'r');
% figure;
% plot(10*log10(error.^2));
err_ensemble_kernel(iter,:) = error .^ 2 ;
system_out_ensemble(iter,:) = system_output .^ 2 ;
model_out_ensemble(iter,:) = model_output .^ 2 ;
multiplier_1_array(iter) = multiplier1 ;
end
plot(mean(system_out_ensemble));
hold on;
plot(mean(model_out_ensemble),'r');
figure;
plot(10 * log10(mean(err_ensemble_kernel)),'r') ;
figure ;

x = (-0.5:0.02:0.5) ;
for i = 1:length(x)
    u = (x(i) / del_x) - floor(x(i) / del_x);
    index = floor(x(i) / del_x)+(Q-1)/2;
    parameter_array_u = [u^2,u,1];
    model_output1 = mean(multiplier_1_array) * parameter_array_u * C * control_point_array1(index:index + 2, 1);
    y(i) = model_output1 ;
    z(i) = f(x(i));
end
plot(x,y) ; hold on ; plot(x,z,'r') ;

function physical_output = f(x)
    k1 = 4 ; k2 = 3 ; h1 = 5 ; h2 = 0.5 ; c1 = -0.8 ; c2 = 0.5 ;
    physical_output = k1 * exp(-1 * ((x - c1) ^ 2) / (2 * h2^2)) +  k2 * exp(-1 * ((x - c2) ^ 2) / (2 * h1^2));
end 

