function [m, k, l, Xo, Yo, Zo, Uo, Vo, Wo]= read_input(input_filename, exp_num)
%The function read_input takes the data from a text file and inputs the
%data into varaibles.
%read_input first determines whether the experiment number is on the list.
%After determining if it is, read_input stores the data (in the same row as
%the experiment number) into several different variables.
A = 1:18;
if ismember(exp_num, A)
classdata = importdata(input_filename);
measurements = classdata.data(exp_num, 2:10);
global m k l Xo Yo Zo Uo Vo Wo;
m = measurements(1);
k = measurements(2);
l = measurements(3);
Xo = measurements(4);
Yo = measurements(5);
Zo = measurements(6);
Uo = measurements(7);
Vo = measurements(8);
Wo = measurements(9);
else
  disp 'Error, invalid input'
m = 'NaN';
k = 'NaN';
l = 'NaN';
Xo = 'NaN';
Yo = 'NaN';
Zo = 'NaN';
Uo = 'NaN';
Vo = 'NaN';
Wo = 'NaN';
end

end