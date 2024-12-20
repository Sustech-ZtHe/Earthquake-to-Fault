function [R,C,scale,minPoi,Mc,hypo]=Detect_Region(file_path)
    hypo =load(file_path);
    if size(hypo,2)>11
        error('Please Check Input File Format!')
    end
    [counts, edges]=histcounts(hypo(:,11));
    [maxCount, maxIndex] = max(counts);
    Mc=edges(maxIndex + 1);minPoi=10;
    lon_dis= max(hypo(:,9))-min(hypo(:,9));
    if lon_dis<=0.1
        R=[0.01,0.1,0.2];C=[1,1.5,2,2.5,3,3.5,4,5];scale=0.1;
    elseif 0.1<lon_dis && lon_dis<=1
        R=[0.1,0.5,1];C=[5,10,15,20,25,30,35,40];scale=1;
    else
        R=[3,6,10];C=[50,60,70,80,90,100,110,120];scale=50;
    end
end
