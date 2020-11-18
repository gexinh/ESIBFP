%import BEM as 'bem'
%import cortex of S1_noise_15db as 'cs'
%% adjust orientation 
chan=cs.GoodChannel;
gainmatrix=bst_gain_orient(bem.Gain(chan,:),bem.GridOrient);
save gainm gainmatrix