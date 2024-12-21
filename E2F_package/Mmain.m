function [a,b,Y,M]=Mmain(ftest4,Mc,hypotout)
    if size(ftest4,1)<10 || size(find(ftest4(:,4)>Mc),1)<10
        M=Mc;
        a=nan;b=nan;Y=nan;
    else
        dt=(max(ftest4(:,4))-min(ftest4(:,4)))/size(ftest4,1);
        Y=[];Y2=[];
        if size(ftest4,1)>500000
            n=1:size(ftest4,1);
            st=size(ftest4,1)/5;
            for i=1:5
                n1=(i-1)*st+1:i*st;
                Y2=[Y2;sum(ftest4(:,4)'>(Mc+dt.*n1)',2)];
            end
            Y=[(Mc+dt.*n)',Y2];
        else
            Y(:,1)=Mc+dt.*(1:size(ftest4,1));
            Y(:,2)=sum(ftest4(:,4)>Mc+dt.*(1:size(ftest4,1)))';
        end
        fun = @(params, x) 10.^(params(1)+params(3).* (params(2) -Y(:,1)));
        params0 = [-1, max(ftest4(:,4)), 1];
        la=[-inf,Mc,0];
        lb=[2,max(hypotout(:,11))+1.2,2];
        % 进行拟合
        params = lsqcurvefit(fun, params0, Y(:,1), Y(:,2),la,lb);
        % 提取拟合参数
        a = params(1);
        M = params(2);
        b = params(3);
    end
end


