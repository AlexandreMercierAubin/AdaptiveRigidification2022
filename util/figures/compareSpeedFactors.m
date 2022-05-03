function [timesFactors,time1,time2] = compareSpeedFactors(td)
    for i = 1:numel(td)/2
        time1(i)=sum(td{2*i-1}.log(2,:));
        time2(i)=sum(td{2*i}.log(2,:));
        timesFactors(i) = time1(i)./time2(i);
    end
end