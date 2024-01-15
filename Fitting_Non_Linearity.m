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
system_output = filter(numerator_coefficients, denominator_coefficients, input);

Q = 21 ; P = 2 ; del_x = 0.2 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;

u = zeros(1,input_length);

index_i1 = ones(1,input_length);
index_i2 = ones(1,input_length);
index_i3 = ones(1,input_length);
index_i4 = ones(1,input_length);
index_i5 = ones(1,input_length);
index_i6 = ones(1,input_length);
index_i7 = ones(1,input_length);
index_i8 = ones(1,input_length);
index_i9 = ones(1,input_length);
index_i10 = ones(1,input_length);

multiplier1 = 1 ;
multiplier2 = 1 ;
multiplier3 = 1 ;
multiplier4 = 1 ;
multiplier5 = 1 ;
multiplier6 = 1 ;
multiplier7 = 1 ;
multiplier8 = 1 ;
multiplier9 = 1 ;
multiplier10 = 1 ;

error = zeros(1,input_length);

model_output = zeros(1,input_length);
model_output1 = zeros(1,input_length);
model_output2 = zeros(1,input_length);
model_output3 = zeros(1,input_length);
model_output4 = zeros(1,input_length);
model_output5 = zeros(1,input_length);
model_output6 = zeros(1,input_length);
model_output7 = zeros(1,input_length);
model_output8 = zeros(1,input_length);
model_output9 = zeros(1,input_length);
model_output10 = zeros(1,input_length);

parameter_array_u = zeros(1,P+1);
parameter_array_u_der = zeros(1,P+1);

control_point_array1 = (-2.5:0.025:-2)';
control_point_array2 = (-2.0:0.025:-1.5)';
control_point_array3 = (-1.5:0.025:-1)';
control_point_array4 = (-1:0.025:-0.5)';
control_point_array5 = (-0.5:0.025:0)';
control_point_array6 = (0:0.025:0.5)';
control_point_array7 = (0.5:0.025:1)';
control_point_array8 = (1:0.025:1.5)';
control_point_array9 = (1.5:0.025:2)';
control_point_array10 = (2:0.025:2.5)';

for i = 1:input_length
    u(i) = (input(i)/del_x) - floor(input(i)/del_x);

    index_i1(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i2(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i3(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i4(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i5(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i6(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i7(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i8(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i9(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i10(i) = floor(input(i)/del_x)+(Q-1)/2;

    parameter_array_u = [(u(i))^2,u(i),1];
    parameter_array_u_der = [2*u(i),1,0];

    model_output1(i) = multiplier1 * parameter_array_u * C * control_point_array1(index_i1(i):index_i1(i)+2,1) ;
    model_output2(i) = multiplier2 * parameter_array_u * C * control_point_array2(index_i2(i):index_i2(i)+2,1) ;
    model_output3(i) = multiplier3 * parameter_array_u * C * control_point_array3(index_i3(i):index_i3(i)+2,1) ;
    model_output4(i) = multiplier4 * parameter_array_u * C * control_point_array4(index_i4(i):index_i4(i)+2,1) ;
    model_output5(i) = multiplier5 * parameter_array_u * C * control_point_array5(index_i5(i):index_i5(i)+2,1) ;
    model_output6(i) = multiplier6 * parameter_array_u * C * control_point_array6(index_i6(i):index_i6(i)+2,1) ;
    model_output7(i) = multiplier7 * parameter_array_u * C * control_point_array7(index_i7(i):index_i7(i)+2,1) ;
    model_output8(i) = multiplier8 * parameter_array_u * C * control_point_array8(index_i8(i):index_i8(i)+2,1) ;
    model_output9(i) = multiplier9 * parameter_array_u * C * control_point_array9(index_i9(i):index_i9(i)+2,1) ;
    model_output10(i) = multiplier10 * parameter_array_u * C * control_point_array10(index_i10(i):index_i10(i)+2,1) ;
    
    model_output(i) = model_output1(i) + model_output2(i) + model_output3(i) + model_output4(i) + model_output5(i) + model_output6(i) + model_output7(i) + model_output8(i) + model_output9(i) + model_output10(i) + system_noise(i) ; 

    error(i) = system_output(i) - model_output(i) ;

    control_point_array1(index_i1(i):index_i1(i)+2 , 1) = control_point_array1(index_i1(i):index_i1(i) + 2 , 1) + 0.025 * multiplier1 * error(i) * C' * parameter_array_u' ;
    control_point_array2(index_i2(i):index_i2(i)+2 , 1) = control_point_array2(index_i2(i):index_i2(i) + 2 , 1) + 0.025 * multiplier2 * error(i) * C' * parameter_array_u' ;
    control_point_array3(index_i3(i):index_i3(i)+2 , 1) = control_point_array3(index_i3(i):index_i3(i) + 2 , 1) + 0.025 * multiplier3 * error(i) * C' * parameter_array_u' ;
    control_point_array4(index_i4(i):index_i4(i)+2 , 1) = control_point_array4(index_i4(i):index_i4(i) + 2 , 1) + 0.025 * multiplier4 * error(i) * C' * parameter_array_u' ;
    control_point_array5(index_i5(i):index_i5(i)+2 , 1) = control_point_array5(index_i5(i):index_i5(i) + 2 , 1) + 0.025 * multiplier5 * error(i) * C' * parameter_array_u' ;
    control_point_array6(index_i6(i):index_i6(i)+2 , 1) = control_point_array6(index_i6(i):index_i6(i) + 2 , 1) + 0.025 * multiplier6 * error(i) * C' * parameter_array_u' ;
    control_point_array7(index_i7(i):index_i7(i)+2 , 1) = control_point_array7(index_i7(i):index_i7(i) + 2 , 1) + 0.025 * multiplier7 * error(i) * C' * parameter_array_u' ;
    control_point_array8(index_i8(i):index_i8(i)+2 , 1) = control_point_array8(index_i8(i):index_i8(i) + 2 , 1) + 0.025 * multiplier8 * error(i) * C' * parameter_array_u' ;
    control_point_array9(index_i9(i):index_i9(i)+2 , 1) = control_point_array9(index_i9(i):index_i9(i) + 2 , 1) + 0.025 * multiplier9 * error(i) * C' * parameter_array_u' ;
    control_point_array10(index_i10(i):index_i10(i)+2 , 1) = control_point_array10(index_i10(i):index_i10(i) + 2 , 1) + 0.025 * multiplier10 * error(i) * C' * parameter_array_u' ;

    multiplier1 = multiplier1 + 0.025 * error(i) * parameter_array_u * C * control_point_array1(index_i1(i):index_i1(i)+2,1) ;
    multiplier2 = multiplier2 + 0.025 * error(i) * parameter_array_u * C * control_point_array2(index_i2(i):index_i2(i)+2,1) ;
    multiplier3 = multiplier3 + 0.025 * error(i) * parameter_array_u * C * control_point_array3(index_i3(i):index_i3(i)+2,1) ;
    multiplier4 = multiplier4 + 0.025 * error(i) * parameter_array_u * C * control_point_array4(index_i4(i):index_i4(i)+2,1) ;
    multiplier5 = multiplier5 + 0.025 * error(i) * parameter_array_u * C * control_point_array5(index_i5(i):index_i5(i)+2,1) ;
    multiplier6 = multiplier6 + 0.025 * error(i) * parameter_array_u * C * control_point_array6(index_i6(i):index_i6(i)+2,1) ;
    multiplier7 = multiplier7 + 0.025 * error(i) * parameter_array_u * C * control_point_array7(index_i7(i):index_i7(i)+2,1) ;
    multiplier8 = multiplier8 + 0.025 * error(i) * parameter_array_u * C * control_point_array8(index_i8(i):index_i8(i)+2,1) ;
    multiplier9 = multiplier9 + 0.025 * error(i) * parameter_array_u * C * control_point_array9(index_i9(i):index_i9(i)+2,1) ;
    multiplier10 = multiplier10 + 0.025 * error(i) * parameter_array_u * C * control_point_array10(index_i10(i):index_i10(i)+2,1) ;
end
err_ensemble_kernel(iter,:) = error .^ 2 ;
system_out_ensemble(iter,:) = system_output .^ 2 ;
model_out_ensemble(iter,:) = model_output .^ 2 ;
end
plot(10 * log10(mean(err_ensemble_kernel)),'b') ; 

x = (-1:0.02:1) ;
for i = 1:length(x)
    u = (x(i) / del_x) - floor(x(i) / del_x);
    index = floor(x(i) / del_x)+(Q-1)/2;
    parameter_array_u = [u^2,u,1];
    model_output1 = multiplier1 * parameter_array_u * C * control_point_array1(index:index + 2, 1);
    model_output2 = multiplier2 * parameter_array_u * C * control_point_array2(index:index + 2, 1);
    model_output3 = multiplier3 * parameter_array_u * C * control_point_array3(index:index + 2, 1);
    model_output4 = multiplier4 * parameter_array_u * C * control_point_array4(index:index + 2, 1);
    model_output5 = multiplier5 * parameter_array_u * C * control_point_array5(index:index + 2, 1);
    model_output6 = multiplier6 * parameter_array_u * C * control_point_array6(index:index + 2, 1);
    model_output7 = multiplier7 * parameter_array_u * C * control_point_array7(index:index + 2, 1);
    model_output8 = multiplier8 * parameter_array_u * C * control_point_array8(index:index + 2, 1);
    model_output9 = multiplier9 * parameter_array_u * C * control_point_array9(index:index + 2, 1);
    model_output10 = multiplier10 * parameter_array_u * C * control_point_array10(index:index + 2, 1);
    y(i) = model_output1 + model_output2 + model_output3 + model_output4 + model_output5 + model_output6 + model_output7 + model_output8 + model_output9 + model_output10; 
end
z = filter(numerator_coefficients, denominator_coefficients, x);
hold on ; 
%plot(x,y,'b') ; hold on ; plot(x,z,'r') ;


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


