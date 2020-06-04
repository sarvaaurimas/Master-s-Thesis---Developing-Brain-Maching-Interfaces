templates = readNPY('templates.npy');
template_1 = templates(75, :, :);
template_1 = reshape(template_1, [82, 64]);
plot(template_1)
