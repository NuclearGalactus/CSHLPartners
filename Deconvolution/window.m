function output = window(in)
    H = hann(length(in));
    output = in .* H;

