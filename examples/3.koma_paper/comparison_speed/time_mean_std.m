clc
KomaCPU = [4.625960;
           4.545896;
           4.555217;
           4.540424;
           4.548038];
KomaGPU0 = [0.198365;
            0.160010;
            0.249766;
            0.179264;
            0.176073];
KomaGPU1 = [0.439456;
            0.349425;
            0.356559;
            0.336995;
            0.354275];
fprintf("Koma CPU = %f +- %f\n", mean(KomaCPU), std(KomaCPU))
fprintf("Koma GPU0 = %f +- %f\n", mean(KomaGPU0), std(KomaGPU0))
fprintf("Koma GPU1 = %f +- %f\n", mean(KomaGPU1), std(KomaGPU1))
%%
MRiLabCPU = [1.525039;
             1.489365;
             1.480434;
             1.506130;
             1.573361];
MRiLabGPU0 = [0.899926;
              0.886749;
              0.961163;
              0.918400;
              0.921536];
MRiLabGPU1 = [0.937016;
              0.868657;
              0.835377;
              0.908583;
              0.838556];
fprintf("MRiLab CPU = %f +- %f\n", mean(MRiLabCPU), std(MRiLabCPU))
fprintf("MRiLab GPU0 = %f +- %f\n", mean(MRiLabGPU0), std(MRiLabGPU0))
fprintf("MRiLab GPU1 = %f +- %f\n", mean(MRiLabGPU1), std(MRiLabGPU1))