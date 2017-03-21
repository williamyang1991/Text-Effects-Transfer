sty = 'flame';
src = 'shu';
trg = 'huo';

imgpath = 'imgs/';

textEffectFinal = text_stylization(sty, src, trg, imgpath);   

imwrite(textEffectFinal, [imgpath, trg, '-', sty, '-', src, '.png']);
