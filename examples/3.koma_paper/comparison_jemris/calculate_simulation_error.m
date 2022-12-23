paths = ["./a.column1d/" "./b.spheresCS/" "./c.brain_suscep/" ...
         "./d.brain_motion/" "./e.brain_spiral/"];
i = 1;
for path = paths
    M_jemris = h5read(path+"signals.h5", sprintf('/signal/channels/%02i',0));
    S_jemris = M_jemris(1,:) + 1i*M_jemris(2,:);
    max_abs_signal = max(abs(S_jemris(:)));
    S_jemris = S_jemris(:);
    
    S_koma   = load(path+'signal_koma.mat').signal;
    
    mae = mean(abs(S_jemris - S_koma)) / max_abs_signal;
    nmrse = norm(S_jemris - S_koma) / norm(S_jemris);

    fprintf("%s\n\tnMAE = %f %%\n\tnRMSE = %f %%\n", path, mae * 100, nmrse * 100)
    
%     figure(i)
%     plot(imag(S_jemris), "b")    
%     hold on
%     plot(imag(S_koma), "r")
%     hold off
%     legend("JEMRIS","Koma")
    
    i = i + 1;
end

figure(i)
plot(imag(S_jemris), "b")    
hold on
plot(imag(S_koma), "r")
plot(abs(S_jemris-S_koma), "k")
hold off
legend("JEMRIS","Koma")
axis([0 1000 -0.02 0.16])