function [] = run_pacing_1AP(BCL, modelswitch)
%% Run a five minute simulation to reach steady-state
sim_length = BCL/1000;
stim = stim_creator_BCL(BCL, 1e-3, sim_length);

if BCL < 1000
    main_human_atrial(stim, ['y0_ss_BCL_0' num2str(BCL) '.dat'], ['1AP_BCL_0' num2str(BCL) '.dat'], sim_length, BCL, 0, modelswitch);
    results_ss_human_atrial(['Results_ss_BCL_0' num2str(BCL)], ['pacing_1AP_BCL_0' num2str(BCL) '.dat']);
else
    main_human_atrial(stim, ['y0_ss_BCL_' num2str(BCL) '.dat'], ['1AP_BCL_' num2str(BCL) '.dat'], sim_length, BCL, 0, modelswitch);
    results_ss_human_atrial(['Results_ss_BCL_' num2str(BCL)], ['pacing_1AP_BCL_' num2str(BCL) '.dat']);
end
