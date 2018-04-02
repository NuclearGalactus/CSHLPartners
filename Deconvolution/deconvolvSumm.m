function output = deconvolvSumm(data, kernelXData, kernel)
halfKern = ceil(length(kernel) / 2) + 1;
dataPadded = [fliplr(data(1:halfKern)) data fliplr(data(end-halfKern+1:end))];
convolved = zeros(length(dataPadded),1);
for i = (halfKern + 1):(length(data))
    currentVal = 0;
    for j = kernelXData
        currentVal = currentVal + dataPadded(i + j) / (kernel(kernelXData == j));
    end
    convolved(i,1) = currentVal;
end
output = convolved;