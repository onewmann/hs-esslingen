function sequence = generateImageSequence(b, v, numFrames)
% Erzeugt eine Bildfolge {b_k(x)} in der angenommen wird:
% b_k(x) = b(x - v*(k - k0)), wobei k0 = ceil(numFrames/2)
[N, M] = size(b);
sequence = zeros(N, M, numFrames);
k0 = ceil(numFrames/2);
for k = 1:numFrames
    shift = v * (k - k0);
    % imtranslate benötigt [dx, dy] = [shift(1), shift(2)]
    frame = imtranslate(b, [shift(1), shift(2)], 'linear', 'FillValues', 0);
    sequence(:, :, k) = frame;
end