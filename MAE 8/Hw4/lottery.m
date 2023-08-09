function [amount] = lottery(ticket)
%Input your lottery ticket to see the amount you've won.
%% The function compares a given ticket to the winning ticket. It then
%% uses the amount of correct numbers to determine the award amount using
%% the function intersect.
load('winning_number.dat');
wticket = winning_number(1:6);
mnumbers = (ticket == wticket);
total = sum(mnumbers, 'all');
if total == 0
    disp('amount = 0');
    amount = 0;
end
if total == 1
	disp('amount = 10');
    amount = 10;
end 
if total == 2
	disp('amount = 100');
    amount = 100;
end 
if total == 3
	disp('amount = 1,000');
    amount = 1000;
end 
if total == 4
	disp('amount = 10,000');
    amount = 10000;
end 
if total == 5
	disp('amount = 1,000,000');
    amount = 1000000;
end
if total == 6
    disp('amount = 100,000,000');
    amount = 100000000;
end
end