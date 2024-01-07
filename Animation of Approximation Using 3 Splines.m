clc;
clear all;
close all;

input_length = 45000 ;
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
index_i2 = ones(1,input_length);
index_i3 = ones(1,input_length);

multiplier1 = 1 ;
multiplier2 = 1 ;
multiplier3 = 1 ;

error = zeros(1,input_length);

model_output = zeros(1,input_length);
model_output1 = zeros(1,input_length);
model_output2 = zeros(1,input_length);
model_output3 = zeros(1,input_length);

parameter_array_u = zeros(1,P+1);
parameter_array_u_der = zeros(1,P+1);

control_point_array1 = (-2:0.2:2)';
control_point_array2 = (-2:0.2:2)';
control_point_array3 = (-2:0.2:2)';

for i = 1:input_length
    i
    u(i) = (input(i)/del_x) - floor(input(i)/del_x);

    index_i1(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i2(i) = floor(input(i)/del_x)+(Q-1)/2;
    index_i3(i) = floor(input(i)/del_x)+(Q-1)/2;

    parameter_array_u = [(u(i))^2,u(i),1];
    parameter_array_u_der = [2*u(i),1,0];

    model_output1(i) = multiplier1 * parameter_array_u * C * control_point_array1(index_i1(i):index_i1(i)+2,1) ;
    model_output2(i) = multiplier2 * parameter_array_u * C * control_point_array2(index_i2(i):index_i2(i)+2,1) ;
    model_output3(i) = multiplier3 * parameter_array_u * C * control_point_array3(index_i3(i):index_i3(i)+2,1) ;
    model_output(i) = model_output1(i) + model_output2(i) + model_output3(i) ; 

    error(i) = system_output(i) - model_output(i) ;

    control_point_array1(index_i1(i):index_i1(i)+2 , 1) = control_point_array1(index_i1(i):index_i1(i) + 2 , 1) + 0.025 * 2 * multiplier1 * error(i) * C' * parameter_array_u' ;
    control_point_array2(index_i2(i):index_i2(i)+2 , 1) = control_point_array2(index_i2(i):index_i2(i) + 2 , 1) + 0.025 * 2 * multiplier2 * error(i) * C' * parameter_array_u' ;
    control_point_array3(index_i3(i):index_i3(i)+2 , 1) = control_point_array3(index_i3(i):index_i3(i) + 2 , 1) + 0.025 * 2 * multiplier3 * error(i) * C' * parameter_array_u' ;
    multiplier1 = multiplier1 + 0.025 * 2 * error(i) * parameter_array_u * C * control_point_array1(index_i1(i):index_i1(i)+2,1) ;
    multiplier2 = multiplier2 + 0.025 * 2 * error(i) * parameter_array_u * C * control_point_array2(index_i2(i):index_i2(i)+2,1) ;
    multiplier3 = multiplier3 + 0.025 * 2 * error(i) * parameter_array_u * C * control_point_array3(index_i3(i):index_i3(i)+2,1) ;
    animate(multiplier1,multiplier2,multiplier3,control_point_array1,control_point_array2,control_point_array3) ;
end

function physical_output = f(x)
    k1 = 4 ; k2 = 3 ; h1 = 5 ; h2 = 0.5 ; c1 = -0.8 ; c2 = 0.5 ;
    physical_output = k1 * exp(-1 * ((x - c1) ^ 2) / (2 * h2^2)) +  k2 * exp(-1 * ((x - c2) ^ 2) / (2 * h1^2));
end 

function graph = animate(multiplier_1,multiplier_2,multiplier3,control_point_array1,control_point_array2,control_point_array3)
    x = (-0.5:0.02:0.5) ;
    Q = 21 ; P = 2 ; del_x = 0.2 ;
    C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;
    fig = gcf;
    for i = 1:length(x)
        u = (x(i) / del_x) - floor(x(i) / del_x);
        index = floor(x(i) / del_x)+(Q-1)/2;
        parameter_array_u = [u^2,u,1];
        model_output1 = multiplier_1 * parameter_array_u * C * control_point_array1(index:index + 2, 1);
        model_output2 = multiplier_2 * parameter_array_u * C * control_point_array2(index:index + 2, 1);
        model_output3 = multiplier3 * parameter_array_u * C * control_point_array3(index:index + 2, 1);
        y(i) = model_output1 + model_output2 + model_output3 ;
        z(i) = f(x(i));  
    end
    clf(fig);
    plot(x,y) ; hold on ; plot(x,z,'r') ;
    pause(1);
end 
