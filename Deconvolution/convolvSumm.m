function output = convolvSumm(data, kernelXData, kernel)
convolved = zeros(length(data),1);
for i = 1:(length(data))
    currentVal = 0;
    for j = kernelXData
        if(i + j > 0 && i + j <= length(data))
        currentVal = currentVal + (kernel(kernelXData == j) * data(i + j));
        end
    end
    convolved(i,1) = currentVal;
end
output = convolved;