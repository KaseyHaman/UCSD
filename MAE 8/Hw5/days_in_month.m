function days = days_in_month(month, leap)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
switch leap
    case 0
        switch month
        case 'jan'
            days = 31;
        case 'feb'
            days = 28;
        case 'mar'
            days = 31;
        case 'apr'
            days = 30;
        case 'may'
            days = 31;
        case 'jun'
            days = 30;
        case 'jul'
            days = 31;
        case 'aug'
            days = 31;
        case 'sep'
            days = 30;
        case 'oct'
            days = 31;
        case 'nov'
            days = 30;
        case 'dec'
            days = 31;
        otherwise 
            days = 'Invalid inputs';
        end 
    case 1
        switch month
        case 'jan'
            days = 31;
        case 'feb'
            days = 29;
        case 'mar'
            days = 31;
        case 'apr'
            days = 30;
        case 'may'
            days = 31;
        case 'jun'
            days = 30;
        case 'jul'
            days = 31;
        case 'aug'
            days = 31;
        case 'sep'
            days = 30;
        case 'oct'
            days = 31;
        case 'nov'
            days = 30;
        case 'dec'
            days = 31;
        otherwise 
            days = 'Invalid inputs';
        end 
end