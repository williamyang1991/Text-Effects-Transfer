startup;

sty = imread('imgs/shu-flame.png');
src = 'shu';
srctext = imread('imgs/shu-text.png');
trgtext = imread('imgs/huo-text.png');

textEffectFinal = text_stylization(sty, src, srctext, trgtext); 

imwrite(textEffectFinal, [imgpath, trg, '-', sty, '-', src, '.png']);

