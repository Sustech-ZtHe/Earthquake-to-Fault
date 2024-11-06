function [R,C,scale,minPoi,Mc,na,hypo]=Detect_Region(file_path)
    hypo =load(file_path);
    %mag
    [counts, edges]=histcounts(hypo(:,11));
    [maxCount, maxIndex] = max(counts);
    Mc=edges(maxIndex + 1);
    na=2;minPoi=10;
    lon_dis= max(hypo(:,9))-min(hypo(:,9));
    if lon_dis<=0.1
        R=[0.01,0.1,0.2];C=[1,2,3,4,5,6,7,8;
                                       1,2,3,4,9,10,11,12;
                                       5,6,7,8,13,14,15,16];scale=0.1;
    elseif 0.1<lon_dis && lon_dis<=1
        R=[0.1,0.5,1];C=[5,8,10,12,15,18,20,25;
                                       1,2,3,4,9,10,11,12;
                                       5,6,7,8,13,14,15,16];scale=1;
    else
        R=[3,6,10];C=[15,20,25,30,35,40,45,50;
                                       1,2,3,4,9,10,11,12;
                                       5,6,7,8,13,14,15,16];scale=50;
    end
end
