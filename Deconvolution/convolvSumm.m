function output = convolvSumm(data2, kernelXData, kernel)
halfKern = ceil(length(kernel) / 2) + 1;
dataPadded = [fliplr(data2(1:halfKern)) data2 fliplr(data2(end-halfKern+1:end))];
convolved = zeros(length(dataPadded),1);
for i = (halfKern + 1):(length(dataPadded) - 2 * halfKern)
    currentVal = 0;
    for j = kernelXData
        currentVal = currentVal + dataPadded(i + j) * (kernel(kernelXData == j));
    end
    convolved(i,1) = currentVal;
end
output = convolved;