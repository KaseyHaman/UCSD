function [total_score, letter] = assign_grade(homework,midterm, project, final)
%This function inputs the scores from the students and outputs a letter
%grade.
%   This function calculates the numerical grade from the student by
%   choosing the function that outputs the highest value. From here, the
%   function sorts the numerical grade into its respective letter grade.
a = 0.25 .* (sum(homework)-min(homework))./7 + 0.25 .* midterm + 0.2 .* project + 0.3 .* final;
b = 0.25 .* (sum(homework)-min(homework))./7 + 0.2 .* project + 0.55 .* final;

if b<a 
    total_score = a;
elseif b>a
    total_score = b;
else
    total_score = a;
end

if total_score >= 93
    letter = 'A';
elseif total_score >= 90
    letter = 'A-';
elseif total_score >= 87
    letter = 'B+';
elseif total_score >= 83
    letter = 'B';
elseif total_score >= 80
    letter = 'B-';
elseif total_score >= 77
    letter = 'C+';
elseif total_score >= 73
    letter = 'C';
elseif total_score >= 70
    letter = 'C-'; 
elseif total_score >= 60
    letter = 'D';
else 
    letter = 'F';

end