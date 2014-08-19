clear stim_1Hz

time_ssajo = 1;
stim_1Hz = stimcreator(1, 1e-3,time_ssajo);

main_human_atrial(stim_1Hz, 'y0_1Hz_vCaNassIk.dat', 'pacing_1AP.dat', 1, 1, 0);

[t, camean] = plot_human_atrial('pacing_1AP.dat');
