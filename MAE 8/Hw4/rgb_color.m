function color = rgb_color(rgb)
%rgb_color takes the input of colors in terms of numbers in order to
%display a range colors
%Call function: rgb_color([rgb]); put a 0 if the color is not present and 
%a 1 if it is: ex) ([red, green, blue]) = ([1, 1, 1]).
if rgb == [1 1 1]
    color = 'white';
elseif rgb == [1 0 1]
    color = 'magenta';
elseif rgb == [1 1 0]
    color = 'yellow';
elseif rgb == [0 1 1]
    color = 'cyan';
elseif rgb == [1 0 0]
    color = 'red';
elseif rgb == [0 1 0]
    color = 'green';
elseif rgb == [0 0 1]
    color = 'blue';
else 
    color = 'Invalid input';
end
end