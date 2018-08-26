startup;

textEffectFinal = text_stylization(imread('imgs/shu-flame.png'), 'shu', imread('imgs/shu-text.png'), imread('imgs/huo-text.png')); 

imwrite(textEffectFinal, [imgpath, trg, '-', sty, '-', src, '.png']);

