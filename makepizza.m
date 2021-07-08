function [inputwithoutsensitive, input, finalsensitive] = makepizza(sample_size, radius1, radius2)
angle1 = 2*pi*rand(sample_size,1);
randomradius1 = radius1*(0.95 + 0.1*rand(sample_size,1));
first = cat(2,randomradius1.*cos(angle1),randomradius1.*sin(angle1),zeros(sample_size,1));

angle2 = 2*pi*rand(sample_size,1);
ramdomradius2 = radius2*(0.95 + 0.1*rand(sample_size,1));
second = cat(2,ramdomradius2.*cos(angle2),ramdomradius2.*sin(angle2),ones(sample_size,1));

concatenated = cat(1,first,second);
idx = randperm(sample_size*2);
input = concatenated(idx,:);
inputwithoutsensitive = input(:,1:2);
finalsensitive = input(:,3);
end 