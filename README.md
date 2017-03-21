# Text-Effects-Transfer

This is a Matlab implementation of the paper.

Shuai Yang, Jiaying Liu, Zhouhui Lian and Zongming Guo, 
Awesome Typography: Statistics-Based Text Effects Transfer, 
Accepted by IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2017.

This code builds upon the PatchMatch algorithm (https://github.com/jbhuang0604/StructCompletion) implemented by Jia-Bin Huang.

It is provided for educational/research purpose only. Please consider citing our paper if you find the software useful for your work.

To run the code, please use the main function text_stylization.m.<br> 
##### Example: 
> textEffectFinal = text_stylization('flame', 'shu', 'huo', 'imgs/', optS); 

   
### External codes

> 1. Flann: for fast approximate nearest neighbor searching.
>>    http://www.cs.ubc.ca/research/flann/

> 2. mirt2D_mexinterp: for fast 2D linear interpolation.
>>   http://www.mathworks.com/matlabcentral/fileexchange/24183-2d-interpolation/content/mirt2D_mexinterp/mirt2D_mexinterp.m


### Contact

Shuai Yang

williamyang@pku.edu.cn
