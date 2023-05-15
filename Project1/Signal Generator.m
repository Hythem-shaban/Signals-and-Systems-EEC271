clear
clc
disp('Welcome To Our Program')
disp('Please ,Enter The Following Specification')
fprintf('\n --------------------------------\n\n')

% request the sampling freq from the user
f=input('Sampling Frequency >> ');

fprintf('\n --------------------------------\n\n')

% request the start and the end of time scale
T_initial = input('Starting Of Time >> ');
T_final = input('Ending Of Time >> ');

fprintf('\n --------------------------------\n\n')

% request the number of breakpints from the user
N=input('Number Of Breakpoints >> ');

Time = [T_initial T_final];

% A loop taking the time at which the breakpoints exist and storing it in
% the vector 'Time'
for i = 1 : 1 : N
    fprintf('The Position of BreakPoint (%d) is at t',i)
    B = input('= ');
    Time = [Time B];
end

% sorting the instants of the breakpoints in the vector 'Time'
Time = sort(Time);
fprintf('\n --------------------------------\n\n')

 
Signal = [];

% These Variables are used in operation part as each list will be useful in
% case of compression and expanion
% A_list store amplitude of all signals (not the slope of ramp)
% M_list store slope of ramp signal
% C_list store intercept of signal
% exponent_list store exponent constant
% freq_list store frequency of sinusodial signal
% theta_list store phase shift of sinusodial signal
function_kind_list = [];
A_list = [];
M_list = [];
C_list = [];
exponent_list = [];
freq_list = [];
theta_list = [];

% A loop for specifying the duration of each function.

% the number of the durations is less than the number of breakpoints by 1
% FOR EX. if the start=0 , end=10, and there are 2 breakpoints at 3, 7
% so there are 4 breakpoints but 3 durations
% from 0 to 3, from 3 to 7, and from 7 to 10 
for i = 1 : 1 : length(Time) - 1 
    Duration = Time(i+1) - Time(i);
    Number_of_Samples = Duration * f;
    t = linspace(Time(i), Time(i+1), Number_of_Samples);
    fprintf('For the interval from %d : %d Choose the function\n',Time(i), Time(i+1))
    
    % menu is a function that displays a modal menu dialog box containing
    % the text in "The_Function" and the choices specified by DC, Ramp,...
    % The menu function returns the number of the selected menu item,
    % or 0 if the user clicks the close button on the window. 
    function_kind = menu('The_Function','DC','Ramp','General Order Polynomial','Exponential','Sinusoidal');
    
    switch  function_kind
        
        % menu returns the number of 'DC' which is 1 so,
        % case 1 is for DC signal which requires the amplitude
        case 1
            disp('You Choosed DC Signal')
            A = input('Enter The Amplitude >> ');
            y = A * ones(1, Number_of_Samples);
            A_list = [A_list A];
        % menu returns the number of 'Ramp' which is 2 so,
        % case 2 is for ramp signal which requires the slope and intercept
        case 2
            disp('You Choosed Ramp Signal')
            M = input('Enter The Slope >> ');
            C = input('Enter The Intercept >> ');
            y = M * t + C;  
            M_list = [M_list M];
            C_list = [C_list C];
            
        % menu returns the number of 'General Order Polynomial' which is 3 so,    
        % case 3 is for General Order Polynomial Signal which requires the
        % amplitude, the power, and the intercept
        case 3
            y = 0;
            disp('You Choosed General Order Polynomial Signal')
            % enter the degree of the polynomial
            power = input('Enter The Power of The Polynomial >> ');
            % create a vector of the amplitudes of each term 
            A = ones(1, power);
            C = input('Enter The Intercept >> ');
            % A loop for assigning the values of the amplitudes and generate the function
            % for ex. if power=2, the func is y = A2*t^2 + A1*t^1 + c
            for j = power : -1 : 1
                fprintf('Enter the Amplitude of t^(%d)', j);
                A(j) = A(j) * input(' >> ');
                y = y + A(j).*(t.^j)  ;
                A_list = [A_list A(j)];
            end
            y = y + C;
            C_list = [C_list C];
        % menu returns the number of 'Exponential' which is 4 so,
        % case 4 is for Exponential Signal which requires the amplitude and
        % the exponent
        case 4
            disp('You Choosed Exponential Signal')
            A = input('Enter The Amplitude >> ');
            Exponent = input('Enter The Exponent >> ');
            y = A * exp(Exponent*t);
            A_list = [A_list A];
            exponent_list = [exponent_list Exponent];
         
        % menu returns the number of 'Sinusoidal' which is 5 so,
        % case 5 is for Sinusoidal Signal which requires the amplitude, the
        % freq, and the phase shift 
        case 5
            disp('You Choosed Sinusoidal Signal')
            A = input('Enter The Amplitude >> ');
            freq = input('Enter The frequency >> ');
            Theta = input('Enter The Phase >> ');
            y = A * sin(2 * pi * freq * t - Theta);
            A_list = [A_list A];
            freq_list = [freq_list freq];
            theta_list = [theta_list Theta];
    end
    % to concatenate all the functions in the full duration of the signal
    Signal = [Signal y];
    fprintf('\n --------------------------------\n\n')
    % function_kind_list will be useful in expansion and compression as it
    % is required to detect signal in order to expand or compress it
    function_kind_list =[function_kind_list function_kind];
 
end
Signal_Samples = (T_final - T_initial) * f;
Signal_Time = linspace(T_initial, T_final, Signal_Samples);
figure;
plot(Signal_Time ,Signal)
% -------------------------------------------------------------------------


% Operation on Signal part 

%load 'test_signal.mat' %this Signal was used for testing Operation on
%Signal Part

% This Part is important as every variable will have effect on the signal
% operation
% Signal_New is used to make operation on this Variable without affecting
% input Signal itself
% flag used as indication. if flag = 0 , program asks user which operation
% would like to perform. if flag = 1, program terminate. this occur only in
% case the user choose None.
% Shift to store Shift value from Shift Operation to be used again
% Factor to store value of compression and expansion value to be used again
% Factor is also used to Store Time Reversing Value (-1)^n , where n is
% number of time Time Reversal has been performed
% Amplitude to store Amplitude Scaling Value

Signal_New = Signal;
flag = 0;
Shift = 0;
Factor = 1;
Amplitude = 1;

while(flag == 0)
disp('Which Operation would you like to use ?')


% this Variables are used as index for lists initialized earlier
L=1;P=1;Q=1;R=1;S=1;U=1;

Operation_kind = menu('Operation on The Signal','Amplitude Scaling','Time Reversal','Time Shift','Expanding the signal','Compressing the signal','None');

% this Part is similar to menu part of the previous part (Signal Generation)
% the first input string is title of the menu
% the next input strings are numbered in ascending order
% Where Amplitude Scaling takes 1 and Time Reversal taks 2 and the rest
% follow the same arrangment

switch Operation_kind
    
    % When menu return number of Amplitude Scaling which is 1
    % switch enter case 1 where amplitude scaling factor will be requested
    
    case 1
        disp('You chose Amplitude Scaling')
        Amplitude=input('Enter Amplitude Scale Factor of new signal >> ');
        Signal_New  = Amplitude.*Signal_New;
        Signal_New_list = [Signal_New];
        
    % When menu return number of Time Reversal which is 2
    % switch enter case 2 where nothing is Requested as signal time will be
    % multiplied by -1 causing reflection of signal in time domain
    
    case 2
        disp('You chose Time Reversal')
        Signal_Time = -1.*Signal_Time;
        Time = -1*Time;
        Factor = -1*Factor;
        Shift = -1*Shift;
        Signal_New_list = [Signal_New];
        
    % When menu return number of Time Shift which is 3
    % switch enter case 3 where Shift in time value is requested.
    % If the input value is positive, the signal is shifted in time to the
    % right.
    % If the input value is negative, the signal is shifted in time to the
    % left. 
    
    case 3
        disp('You chose Time Shift')
        Shift_input = input('Amount of shift in Time (+ve -> Right , -ve <- Left) >> ');
        Time = Time + Shift_input;
        Shift = Shift + Shift_input;
        Signal_Time = linspace(Time(1), Time(end), Signal_Samples);
        Signal_New_list = [Signal_New];
        
    % When menu return number of Time Shift which is 4
    % switch enter case 4 where Expanding Factor is requested.
    % the input value must be greater than 1.
    % the signal will be expanded in time by this factor
    
    case 4
        disp('You chose Expanding the Signal')
        Expand = input('Expanding Factor (Value > 1) >> ');
        while (Expand < 1)
            disp('Expanding Factor must be greater than 1')
            Expand = input('Expanding Factor (Value > 1) >> ');
        end    
        
        % This part is used to identify the type of previously chosen
        % function to apply expansion on it so that number of samples in
        % indpendent Variable Signal_Time equal to that of dependent 
        % Variable Signal_New
        
        Signal_New_list =[];
        Signal_Samples = 0;
        Factor = Factor*Expand;
        Time = Time* Expand;
        Shift = Shift *  Expand;
        for k = 1 : 1 : N+1
            
         % abs is used as Duration must be +Ve if Time Reversal was
         % previously performed the Duration would be -Ve causing error in
         % linspace as Number of Samples should be -Ve
         
         Duration = abs(Time(k+1) - Time(k));
         Number_of_Samples = round(Duration * f);
         Signal_Samples = Number_of_Samples + Signal_Samples;
         Signal_Time = linspace(Time(k), Time(k+1), Number_of_Samples);
         
         % This part is used to redraw the signal by detecting each part
         % and apply changes due to previous input
         
        switch function_kind_list(k)
            case 1
              Signal_New = A_list(L)*ones(1,Number_of_Samples);
              Signal_New_list = [Signal_New_list Signal_New];
              L=L+1;
            case 2
                Signal_New =  M_list(P) *(1/ Factor)*(Signal_Time -Shift) + C_list(Q);
                Signal_New_list = [Signal_New_list Signal_New];
                P=P+1;
                Q=Q+1;
            case 3
                y = 0;
            for j = power : -1 : 1
                y = y + A_list(power-j+L).*(((1/ Factor)*(Signal_Time -Shift)).^j);
            end
            Signal_New = y + C_list(Q);
            Signal_New_list = [Signal_New_list Signal_New];
            L=L+power-1+1;
            Q=Q+1;
            case 4
            Signal_New = A_list(L) * exp(exponent_list(R)*((1/ Factor)*(Signal_Time -Shift)));
            Signal_New_list = [Signal_New_list Signal_New];
            L=L+1;
            R=R+1;
            case 5
            Signal_New = A_list(L) * sin(2 * pi * freq_list(S) *((1/ Factor)*(Signal_Time -Shift)) - theta_list(U)); 
            Signal_New_list = [Signal_New_list Signal_New];
            L=L+1;
            S=S+1;
            U=U+1;
        end
        end
        Signal_Time = linspace(Time(1), Time(end), Signal_Samples);
        Signal_New_list = Signal_New_list*Amplitude;
        
    % When menu return number of Time Shift which is 5
    % switch enter case 5 where Comressing Factor  is requested.
    % the input value must be less than 1.
    % the signal will be compressed in time by this factor
    
    case 5
        disp('You chose Compressing the Signal')
        compress = input('Compressing Factor (Value < 1) >> ');
        while (compress > 1 | compress <= 0)
            disp('Compressing Falue must be less than 1 and greater than 0')
            compress = input('Compressing Factor (Value < 1) >> ');
        end 
        
        % This part is used to identify the type of previously chosen
        % function to apply expansion on it so that number of samples in
        % indpendent Variable Signal_Time equal to that of dependent 
        % Variable Signal_New
        
        Signal_New_list =[];
        Signal_Samples = 0;
        Factor = Factor*compress;
        Time = Time*compress;
        Shift = Shift * compress;
        for k = 1 : 1 : N+1
            
         % abs is used as Duration must be +Ve if Time Reversal was
         % previously performed the Duration would be -Ve causing error in
         % linspace as Number of Samples should be -Ve
         
         Duration = abs(Time(k+1) - Time(k));
         Number_of_Samples = round(Duration * f);
         Signal_Samples = Number_of_Samples + Signal_Samples;
         Signal_Time = linspace(Time(k), Time(k+1), Number_of_Samples);
         
         % This part is used to redraw the signal by detecting each part
         % and apply changes due to previous input
         
        switch function_kind_list(k)
            case 1
              Signal_New = A_list(L)*ones(1,Number_of_Samples);
              Signal_New_list = [Signal_New_list Signal_New];
              L=L+1;
            case 2
                Signal_New =  M_list(P) *(1/Factor)*(Signal_Time -Shift) + C_list(Q);
                Signal_New_list = [Signal_New_list Signal_New];
                P=P+1;
                Q=Q+1;
            case 3
                y = 0;
            for j = power : -1 : 1
                y = y + A_list(power-j+L).*(((1/Factor)*(Signal_Time -Shift)).^j)  ;
            end
            Signal_New = y + C_list(Q);
            Signal_New_list = [Signal_New_list Signal_New];
            L=L+power-1+1;
            Q=Q+1;
            case 4
            Signal_New = A_list(L) * exp(exponent_list(R)*((1/Factor)*(Signal_Time -Shift)));
            Signal_New_list = [Signal_New_list Signal_New];
            L=L+1;
            R=R+1;
            case 5
            Signal_New = A_list(L) * sin(2 * pi * freq_list(S) *((1/Factor)*(Signal_Time -Shift)) - theta_list(U)); 
            Signal_New_list = [Signal_New_list Signal_New];
            L=L+1;
            S=S+1;
            U=U+1;
        end
        end
        Signal_Time = linspace(Time(1), Time(end), Signal_Samples);
        Signal_New_list = Signal_New_list*Amplitude;
        
    % When menu return number of Time Shift which is 5
    % Noting will happen in signal
    
    case 6
        disp('You chose None')
        flag = 1;
     Signal_New_list = [Signal_New];
end
fprintf('\n --------------------------------\n\n')

disp('Signal After Operation')

fprintf('\n --------------------------------\n\n')

% Signal_New must be equal to Signal_New_list to prepare for next operation

Signal_New = [Signal_New_list];

% Displaying the new signal

figure;
plot(Signal_Time,Signal_New_list)

end