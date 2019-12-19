%% Subtract Background
function I = subtract_background(I, bg)
    sel = I <= bg;
    I(sel) = 0;
end