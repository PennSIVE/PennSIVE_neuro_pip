# Required packages:
suppressMessages(library(argparser))
suppressMessages(library(neurobase))
suppressMessages(library(tidyverse))
suppressMessages(library(oro.nifti))
suppressMessages(library(oro.dicom))
suppressMessages(library(WhiteStripe))
suppressMessages(library(fslr))
suppressMessages(library(ANTsR))
suppressMessages(library(extrantsr))
suppressMessages(library(mimosa))
suppressMessages(library(parallel))
suppressMessages(library(purrr))
suppressMessages(library(pbmcapply))
suppressMessages(library(pbapply))
suppressMessages(library(ANTsRCore))
suppressMessages(library(stringr))
suppressMessages(library(caret))


p <- arg_parser("Running Paramagnetic Rim Lesion (PRL) detection pipeline to obtain PRL probability.", hide.opts = FALSE)
p <- add_argument(p, "--mainpath", short = '-m', help = "Specify the main path where MRI images can be found.")
p <- add_argument(p, "--participant", short = '-p', help = "Specify the subject id.")
p <- add_argument(p, "--session", short = '-s', help = "Specify the session id.")
p <- add_argument(p, "--t1", help = "Specify the T1 sequence name.")
p <- add_argument(p, "--t2", help = "Specify the T2 sequence name.")
p <- add_argument(p, "--flair", help = "Specify the FLAIR sequence name.")
p <- add_argument(p, "--phase", help = "Specify the Phase sequence name.")
p <- add_argument(p, "--n4", help = "Specify whether to run bias correction step.", default = TRUE)
p <- add_argument(p, "--skullstripping", short = '-s', help = "Specify whether to run skull stripping step.", default = FALSE)
p <- add_argument(p, "--registration", short = '-r', help = "Specify whether to run registration step.", default = TRUE)
p <- add_argument(p, "--whitestripe", short = '-w', help = "Specify whether to run whitestripe step.", default = TRUE)
p <- add_argument(p, "--mimosa", help = "Specify whether to run mimosa segmentation step.", default = TRUE)
p <- add_argument(p, "--threshold", help = "Specify the threshold used to generate mimosa mask.", default = 0.2)
p <- add_argument(p, "--dilation", short = '-d', help = "Specify whether to dilate lesion.", default = TRUE)
p <- add_argument(p, "--step", help = "Specify the step of PRL pipeline. preparation, PRL_run or consolidation.", default = "preparation")
p <- add_argument(p, "--lesioncenter", help = "Provide the path to the lesioncenter package.")
p <- add_argument(p, "--mpath", help = "Specify the path to the trained mimosa model.")
p <- add_argument(p, "--aprlpath", help = "Specify the path to the trained aprl model.")
p <- add_argument(p, "--helpfunc", help = "Specify the path to the help functions.")
argv <- parse_args(p)

# Read in Files
main_path = argv$mainpath
if (argv$step == "preparation"){
  # Load lesion center package
  my_path = paste0(argv$lesioncenter, "/")
  source_files = list.files(my_path)
  purrr::map(paste0(my_path, source_files), source)
  p = argv$participant
  ses = argv$session
  message('Checking inputs...')
  if(is.na(argv$t1)) stop("Missing T1 sequence!")else{
    t1 = readnii(paste0(main_path, "/data/", p, "/", ses, "/anat/", argv$t1))
  }
  if(is.na(argv$flair)) stop("Missing FLAIR sequence!")else{
    flair = readnii(paste0(main_path, "/data/", p, "/", ses, "/anat/", argv$flair))
  }
  if(is.na(argv$phase)) stop("Missing Phase sequence!")else{
    phase = read_rpi(paste0(main_path, "/data/", p, "/", ses, "/anat/", argv$phase))
  }
  if(!is.na(argv$t2)){
    t2 = readnii(paste0(main_path, "/data/", p, "/", ses, "/anat/", argv$t2))
  }

  # Bias Correction
  bias.out.dir = paste0(main_path, "/data/", p, "/", ses, "/bias_correction")
  if(argv$n4){
    # only create directory if it doesn't already exist - EAH 5/7/26
    if (!dir.exists(bias.out.dir)) {
      dir.create(bias.out.dir,showWarnings = FALSE)
    }
    
    # file check - EAH 5/7/26
    t1_biascorrect_path = paste0(bias.out.dir,"/T1_n4.nii.gz")
    flair_biascorrect_path = paste0(bias.out.dir,"/FLAIR_n4.nii.gz")
    phase_biascorrect_path = paste0(bias.out.dir,"/PHASE_n4.nii.gz")
    t2_biascorrect_path = paste0(bias.out.dir,"/T2_n4.nii.gz")
    
    if (is.na(argv$t2)) {
      files_n4 = c(t1_biascorrect = t1_biascorrect_path, flair_biascorrect = flair_biascorrect_path, phase_biascorrect = phase_biascorrect_path)
      missing_files_n4 = files_n4[!file.exists(c(t1_biascorrect_path, flair_biascorrect_path, phase_biascorrect_path))]
    } else {
      files_n4 = c(t1_biascorrect = t1_biascorrect_path, flair_biascorrect = flair_biascorrect_path, phase_biascorrect = phase_biascorrect_path, t2_biascorrect = t2_biascorrect_path)
      missing_files_n4 = files_n4[!file.exists(c(t1_biascorrect_path, flair_biascorrect_path, phase_biascorrect_path, t2_biascorrect_path))]
    }
    
    if ('flair_biascorrect' %in% names(missing_files_n4)) {
      flair_biascorrect = bias_correct(file = flair,
                                      correction = "N4",
                                      verbose = TRUE)
      writenii(flair_biascorrect,paste0(bias.out.dir,"/FLAIR_n4.nii.gz"))
    } else {
      flair_biascorrect = readnii(flair_biascorrect_path)
    }
    
    if ('t1_biascorrect' %in% names(missing_files_n4)) {
      t1_biascorrect = bias_correct(file = t1,
                                   correction = "N4",
                                   verbose = TRUE)
      writenii(t1_biascorrect,paste0(bias.out.dir,"/T1_n4.nii.gz"))
    } else {
      t1_biascorrect = readnii(t1_biascorrect_path)
    }
    
    if ('phase_biascorrect' %in% names(missing_files_n4)) {
      phase_biascorrect = bias_correct(file = phase,
                                   correction = "N4",
                                   verbose = TRUE)
      writenii(phase_biascorrect,paste0(bias.out.dir,"/PHASE_n4.nii.gz"))
    } else {
      phase_biascorrect = readnii(phase_biascorrect_path)
    }
    
    if(!is.na(argv$t2)){
      if ('t2_biascorrect' %in% names(missing_files_n4)) {
        t2_biascorrect = bias_correct(file = t2,
                                   correction = "N4",
                                   verbose = TRUE)
        writenii(t2_biascorrect,paste0(bias.out.dir,"/T2_n4.nii.gz"))
      } else {
        t2_biascorrect = readnii(t2_biascorrect_path)
      }
    }
  }else{
    t1_biascorrect = readnii(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/T1_n4.nii.gz"))
    flair_biascorrect = readnii(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/FLAIR_n4.nii.gz"))
    if(!is.na(argv$t2)){
      t2_biascorrect = readnii(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/T2_n4.nii.gz")) 
    }
    existing_files = list.files(paste0(main_path, "/data/", p, "/", ses, "/bias_correction"))
    phase_file = existing_files[which(grepl("PHASE_n4*", existing_files))]
    if(length(phase_file) == 0){
      phase_biascorrect = bias_correct(file = phase,
                                 correction = "N4",
                                 verbose = TRUE)
    writenii(phase_biascorrect,paste0(bias.out.dir,"/PHASE_n4.nii.gz"))
    }else{
      phase_biascorrect = read_rpi(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/PHASE_n4.nii.gz"))
    }
  }

  # Skull Stripping
  brain.out.dir = paste0(main_path, "/data/", p, "/", ses, "/t1_brain")
  if(!argv$skullstripping){
    brain_paths = list.files(brain.out.dir, recursive = TRUE, full.names = TRUE)
    brain_path = brain_paths[which(grepl("*brain.nii.gz$", brain_paths))]
    brain_mask_path = brain_paths[which(grepl("*brainmask.nii.gz$", brain_paths))]
    t1_fslbet_robust = readnii(brain_path)
    brain_mask = readnii(brain_mask_path)
    existing_files = list.files(bias.out.dir)
    t1_n4_brain_file = existing_files[which(grepl("T1_brain_n4*", existing_files))]
    if(length(t1_n4_brain_file) == 0){
      t1_fslbet_robust = bias_correct(file = t1_fslbet_robust,
                               correction = "N4",
                               verbose = TRUE)
      # only write file if it doesn't exist - EAH 5/7/26
      if (!file.exists(paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))) {
        writenii(t1_fslbet_robust,paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
      }
    }else{
      t1_fslbet_robust = readnii(paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))
      } 
    }

  if (argv$skullstripping){
    # only create directory if it doesn't already exist - EAH 5/6/26
    if (!dir.exists(brain.out.dir)) {
      dir.create(brain.out.dir,showWarnings = FALSE)
    }
    
    # file check - EAH 5/6/26
    t1_fslbet_robust_path = paste0(brain.out.dir,"/T1_brain.nii.gz")
    brain_mask_path = paste0(brain.out.dir,"/T1_brainmask.nii.gz")
    
    files_skullstripping = c(t1_fslbet_robust = t1_fslbet_robust_path, brain_mask = brain_mask_path)
    missing_files_skullstripping = files_skullstripping[!file.exists(c(t1_fslbet_robust_path, brain_mask_path))]
    
    if ('t1_fslbet_robust' %in% names(missing_files_skullstripping)) {
      t1_fslbet_robust = fslbet_robust(t1_biascorrect,reorient = FALSE,correct = FALSE)
      writenii(t1_fslbet_robust,paste0(brain.out.dir,"/T1_brain.nii.gz"))
      writenii(t1_fslbet_robust,paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
    } else {
      t1_fslbet_robust = readnii(t1_fslbet_robust_path)
    }
    if ('brain_mask' %in% names(missing_files_skullstripping)) {
      brain_mask = t1_fslbet_robust > 0 
      writenii(brain_mask,paste0(brain.out.dir,"/T1_brainmask.nii.gz"))
    } else {
      brain_mask = readnii(brain_mask_path)
    }
  }

  # Registration to FLAIR Space
  reg.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/FLAIR_space")
  if (argv$registration){
    # only create directory if it doesn't already exist - EAH 5/7/26
    if (!dir.exists(reg.out.dir)) {
      dir.create(reg.out.dir,showWarnings = FALSE, recursive = TRUE)
    }
    
    # file check - EAH 5/7/26
    t1_to_flair_path = paste0(reg.out.dir,"/t1_reg_to_flair0GenericAffine.mat")
    t1_reg_path = paste0(reg.out.dir,"/t1_n4_brain_reg_flair.nii.gz")
    brainmask_reg_path = paste0(reg.out.dir,"/brainmask_reg_flair.nii.gz")
    flair_n4_brain_path = paste0(reg.out.dir,"/flair_n4_brain.nii.gz")
    t2_to_flair_path = paste0(reg.out.dir,"/t2_reg_to_flair0GenericAffine.mat")
    t2_n4_brain_path = paste0(reg.out.dir,"/t2_n4_brain_reg_flair.nii.gz")
    
    if (is.na(argv$t2)) {
      files_reg_flair = c(t1_to_flair = t1_to_flair_path, t1_reg = t1_reg_path, brainmask_reg = brainmask_reg_path, flair_n4_brain = flair_n4_brain_path)
      missing_files_reg_flair = files_reg_flair[!file.exists(c(t1_to_flair_path, t1_reg_path, brainmask_reg_path, flair_n4_brain_path))]
    } else {
      files_reg_flair = c(t1_to_flair = t1_to_flair_path, t1_reg = t1_reg_path, brainmask_reg = brainmask_reg_path, flair_n4_brain = flair_n4_brain_path, t2_to_flair = t2_to_flair_path, t2_n4_brain = t2_n4_brain_path)
      missing_files_reg_flair = files_reg_flair[!file.exists(c(t1_to_flair_path, t1_reg_path, brainmask_reg_path, flair_n4_brain_path, t2_to_flair_path, t2_n4_brain_path))]
    }
    
    ## Register T1 to FLAIR space 
    if ('t1_to_flair' %in% names(missing_files_reg_flair)) {
      t1_to_flair = registration(filename = t1_biascorrect,
                               template.file = flair_biascorrect,
                               typeofTransform = "Rigid", remove.warp = FALSE,
                               outprefix=paste0(reg.out.dir,"/t1_reg_to_flair")) 
    } else {
      t1_to_flair = readAntsrTransform(paste0(reg.out.dir,"/t1_reg_to_flair0GenericAffine.mat"))
    }

    if ('t1_reg' %in% names(missing_files_reg_flair)) {
      t1_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(t1_fslbet_robust),
                                          transformlist = t1_to_flair$fwdtransforms, interpolator = "welchWindowedSinc"))
      writenii(t1_reg, paste0(reg.out.dir,"/t1_n4_brain_reg_flair"))}
    else {
      t1_reg = readnii(t1_reg_path)
    }
    
    if ('brainmask_reg' %in% names(missing_files_reg_flair)) {
      brainmask_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(brain_mask),
                                               transformlist = t1_to_flair$fwdtransforms, interpolator = "nearestNeighbor"))
      writenii(brainmask_reg, paste0(reg.out.dir,"/brainmask_reg_flair"))
    } else {
      brainmask_reg = readnii(brainmask_reg_path)
    }
    
    if ('flair_n4_brain' %in% names(missing_files_reg_flair)) {
      flair_n4_brain = flair_biascorrect
      flair_n4_brain[brainmask_reg==0] = 0
      writenii(flair_n4_brain, paste0(reg.out.dir,"/flair_n4_brain"))
    } else {
      flair_n4_brain = readnii(flair_n4_brain_path)
    }

    ## Register T2 to FLAIR space 
    if(!is.na(argv$t2)){
      if ('t2_to_flair' %in% names(missing_files_reg_flair)) {
        t2_to_flair = registration(filename = t2_biascorrect,
                                   template.file = flair_biascorrect,
                                   typeofTransform = "Rigid", remove.warp = FALSE,
                                   outprefix=paste0(reg.out.dir,"/t2_reg_to_flair"))
      } else {
        t2_to_flair = readAntsrTransform(paste0(reg.out.dir,"/t2_reg_to_flair0GenericAffine.mat"))
      }
      if ('t2_n4_brain' %in% names(missing_files_reg_flair)) {
        t2_n4_brain = t2_to_flair$outfile
        t2_n4_brain[brainmask_reg==0] = 0
        writenii(t2_n4_brain, paste0(reg.out.dir,"/t2_n4_brain_reg_flair"))
      } else {
        t2_n4_brain = readnii(t2_n4_brain_path)
      }
    }
  }else{
    t1_reg = readnii(paste0(reg.out.dir, "/t1_n4_brain_reg_flair.nii.gz"))
    flair_n4_brain = readnii(paste0(reg.out.dir, "/flair_n4_brain.nii.gz"))
    brainmask_reg = readnii(paste0(reg.out.dir, "/brainmask_reg_flair"))
    if(!is.na(argv$t2)){
      t2_reg = readnii(paste0(reg.out.dir, "/t2_n4_brain_reg_flair.nii.gz"))
    }
  }

  # WhiteStripe normalize data
  white.out.dir = paste0(main_path, "/data/", p, "/", ses, "/whitestripe/FLAIR_space")
  if(argv$whitestripe){
    # only create directory if it doesn't already exist - EAH 5/7/26
    if (!dir.exists(white.out.dir)) {
      dir.create(white.out.dir,showWarnings = FALSE, recursive = TRUE)
    }
    
    # file check - EAH 5/7/26
    t1_n4_reg_brain_ws_path = paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws.nii.gz")
    t2_n4_reg_brain_ws_path = paste0(white.out.dir, "/t2_n4_brain_reg_flair_ws.nii.gz")
    flair_n4_brain_ws_path = paste0(white.out.dir, "/flair_n4_brain_ws.nii.gz")
    
    if (is.na(argv$t2)) {
      files_ws_flair = c(t1_n4_reg_brain_ws = t1_n4_reg_brain_ws_path, flair_n4_brain_ws = flair_n4_brain_ws_path)
      missing_files_ws_flair = files_ws_flair[!file.exists(c(t1_n4_reg_brain_ws_path, flair_n4_brain_ws_path))]
    } else {
      files_ws_flair = c(t1_n4_reg_brain_ws = t1_n4_reg_brain_ws_path, t2_n4_reg_brain_ws = t2_n4_reg_brain_ws_path, flair_n4_brain_ws = flair_n4_brain_ws_path)
      missing_files_ws_flair = files_ws_flair[!file.exists(c(t1_n4_reg_brain_ws_path, t2_n4_reg_brain_ws_path, flair_n4_brain_path, flair_n4_brain_ws_path))]
    }
    
    if ('t1_n4_reg_brain_ws' %in% names(missing_files_ws_flair)) {
      ind1 = whitestripe(t1_reg, "T1")
      t1_n4_reg_brain_ws = whitestripe_norm(t1_reg, ind1$whitestripe.ind)
      writenii(t1_n4_reg_brain_ws, paste0(white.out.dir,"/t1_n4_brain_reg_flair_ws"))
    } else {
      t1_n4_reg_brain_ws = readnii(t1_n4_reg_brain_ws_path)
    }
    
    
    if(!is.na(argv$t2)){
      if ('t2_n4_reg_brain_ws' %in% names(missing_files_ws_flair)) {
        ind2 = whitestripe(t2_n4_brain, "T2")
        t2_n4_reg_brain_ws = whitestripe_norm(t2_n4_brain, ind2$whitestripe.ind)
        writenii(t2_n4_reg_brain_ws, paste0(white.out.dir,"/t2_n4_brain_reg_flair_ws"))
      } else {
        t2_n4_reg_brain_ws = readnii(t2_n4_reg_brain_ws_path)
      }
    }
    
    if ('flair_n4_brain_ws' %in% names(missing_files_ws_flair)) {
      ind3 = whitestripe(flair_n4_brain, "T2")
      flair_n4_brain_ws = whitestripe_norm(flair_n4_brain, ind3$whitestripe.ind)
      writenii(flair_n4_brain_ws, paste0(white.out.dir,"/flair_n4_brain_ws"))
    } else {
      flair_n4_brain_ws = readnii(flair_n4_brain_ws_path)
    }
    }else{
      t1_n4_reg_brain_ws = readnii(paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws"))
      if(!is.na(argv$t2)){
        t2_n4_reg_brain_ws = readnii(paste0(white.out.dir, "/t2_n4_brain_reg_flair_ws"))
      }
      flair_n4_brain_ws = readnii(paste0(white.out.dir, "/flair_n4_brain_ws"))
    }

  # Mimosa
  mim.out.dir = paste0(main_path, "/data/", p, "/", ses, "/mimosa")
  if(argv$mimosa){
    # only create directory if it doesn't already exist - EAH 5/7/26
    if (!dir.exists(mim.out.dir)) {
      dir.create(mim.out.dir,showWarnings = FALSE)
    }
    
    # file check - EAH 5/7/26
    probmap_path = paste0(mim.out.dir,"/mimosa_prob.nii.gz")
    mimosa_mask_path = paste0(mim.out.dir,"/mimosa_mask.nii.gz")
    
    files_mimosa = c(probmap = probmap_path, mimosa_mask = mimosa_mask_path)
    missing_files_mimosa = files_mimosa[!file.exists(c(probmap_path, mimosa_mask_path))]

    if ('probmap' %in% names(missing_files_mimosa)) {
      mimosa = mimosa_data(brain_mask=brainmask_reg, FLAIR=flair_n4_brain_ws, T1=t1_n4_reg_brain_ws, gold_standard=NULL, normalize="no", cores = 1, verbose = TRUE)
      mimosa_df = mimosa$mimosa_dataframe
      cand_voxels = mimosa$top_voxels
      tissue_mask = mimosa$tissue_mask
      load(argv$mpath) 
      predictions_WS = predict(mimosa_model, mimosa_df, type="response")
      predictions_nifti_WS = niftiarr(cand_voxels, 0)
      predictions_nifti_WS[cand_voxels==1] = predictions_WS
      probmap = fslsmooth(predictions_nifti_WS, sigma = 1.25, mask=tissue_mask, retimg=TRUE, smooth_mask=TRUE) 
      writenii(probmap, paste0(mim.out.dir,"/mimosa_prob"))
    } else {
      probmap = readnii(probmap_path)
    }
    if ('mimosa_mask' %in% names(missing_files_mimosa)) {
      writenii(probmap > as.numeric(argv$threshold), paste0(mim.out.dir,"/mimosa_mask"))
    } 
    
  }else{
    probmap = readnii(paste0(mim.out.dir,"/mimosa_prob"))
  }

  # Register to EPI Space
  reg.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/EPI_space")
  #if(!file.exists(reg.epi.out.dir)){
    # only create directory if it doesn't already exist - EAH 5/7/26
    if (!dir.exists(reg.epi.out.dir)) {
      dir.create(reg.epi.out.dir,showWarnings = FALSE)
    }
    
    # file check - EAH 5/7/26
    
    brainmask_reg_epi_path = paste0(reg.epi.out.dir,'/brainmask_reg_epi.nii.gz')
    phase_n4_brain_path = paste0(reg.epi.out.dir,'/phase_n4_brain.nii.gz')
    t1_reg_epi_path = paste0(reg.epi.out.dir, "/t1_reg_epi.nii.gz")
    flair_reg_epi_path = paste0(reg.epi.out.dir, "/flair_reg_epi.nii.gz")
    mimosa_reg_epi_path = paste0(reg.epi.out.dir, "/mimosa_reg_epi.nii.gz")
    mimosa_mask_reg_epi_path = paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi.nii.gz")
    fast_reg_epi_path = paste0(reg.epi.out.dir, "/fast_reg_epi.nii.gz")
    
    files_reg_epi = c(brainmask_reg_epi = brainmask_reg_epi_path, phase_n4_brain = phase_n4_brain_path, t1_reg_epi = t1_reg_epi_path, flair_reg_epi = flair_reg_epi_path,
                      mimosa_reg_epi = mimosa_reg_epi_path, mimosa_mask_reg_epi = mimosa_mask_reg_epi_path, fast_reg_epi = fast_reg_epi_path)
    missing_files_reg_epi = files_reg_epi[!file.exists(c(brainmask_reg_epi_path, phase_n4_brain_path, t1_reg_epi_path, flair_reg_epi_path,
                                                         mimosa_reg_epi_path, mimosa_mask_reg_epi_path, fast_reg_epi_path))]
    
    # Fast
    t1_fast = fast(t1_reg, bias_correct = F, opts = "-t 1 -n 3")
    flair_to_epi = registration(filename = flair_biascorrect,
                                  template.file = abs(phase_biascorrect),
                                  typeofTransform = "Rigid", remove.warp = FALSE) ### rigid

    if ('brainmask_reg_epi' %in% names(missing_files_reg_epi)) {
      brainmask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_biascorrect)), moving = oro2ants(brainmask_reg),
                                                      transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      writenii(brainmask_reg_epi, paste0(reg.epi.out.dir,'/brainmask_reg_epi'))
    } else {
      brainmask_reg_epi = readnii(brainmask_reg_epi_path)
    }
    
    if ('phase_n4_brain' %in% names(missing_files_reg_epi)) {
      phase_n4_brain = phase_biascorrect * brainmask_reg_epi
      writenii(phase_n4_brain, paste0(reg.epi.out.dir,'/phase_n4_brain'))
    } else {
      phase_n4_brain = readnii(phase_n4_brain_path)
    }

    if ('t1_reg_epi' %in% names(missing_files_reg_epi)) {
      t1_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(t1_reg),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
      writenii(t1_reg_epi, paste0(reg.epi.out.dir, "/t1_reg_epi"))
    } else {
      t1_reg_epi = readnii(t1_reg_epi_path)
    }
    
    if ('flair_reg_epi' %in% names(missing_files_reg_epi)) {
      flair_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(flair_n4_brain),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
      writenii(flair_reg_epi, paste0(reg.epi.out.dir, "/flair_reg_epi"))
    } else {
      flair_reg_epi = readnii(flair_reg_epi_path)
    }
    
    if ('mimosa_reg_epi' %in% names(missing_files_reg_epi)) {
      mimosa_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(probmap),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
      writenii(mimosa_reg_epi, paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
    } else {
      mimosa_reg_epi = readnii(mimosa_reg_epi_path)
    }
    
    if ('mimosa_mask_reg_epi' %in% names(missing_files_reg_epi)) {
      mimosa_mask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(probmap>0.2),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      writenii(mimosa_mask_reg_epi, paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
    } else {
      mimosa_mask_reg_epi = readnii(mimosa_mask_reg_epi_path)
    }
    
    if ('fast_reg_epi' %in% names(missing_files_reg_epi)) {
      fast_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(t1_fast),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      writenii(fast_reg_epi, paste0(reg.epi.out.dir,'/fast_reg_epi'))
    } else {
      fast_reg_epi = readnii(fast_reg_epi_path)
    }
  # }else{
  #   mimosa_reg_epi = readnii(paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
  #   mimosa_mask_reg_epi = readnii(paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
  #   fast_reg_epi = readnii(paste0(reg.epi.out.dir,'/fast_reg_epi'))
  # }

  # Dilating Lesions
  if (argv$dilation){
    if(!file.exists(paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi_dilated.nii.gz"))){
      lesmask_dil = fsldilate(mimosa_mask_reg_epi) # dilate segmentation mask by 1 voxel
      dil = lesmask_dil; dil[mimosa_mask_reg_epi==1] = 0 # get just the dilated voxels
      dil[fast_reg_epi==3] = 0 # find dilated voxels in gm/csf (subset out wm voxels)
      lesmask_dil_mask = lesmask_dil - dil # take out gm/csf dilated voxels
      writenii(lesmask_dil_mask, paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi_dilated"))
    }else{lesmask_dil_mask = readnii(paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi_dilated"))}
  }else{
      lesmask_dil_mask = mimosa_mask_reg_epi
  }

  # Lesion Labeling
  if(!file.exists(paste0(reg.epi.out.dir, "/lesions_reg_epi_labeled.nii.gz"))){
    source(paste0(argv$helpfunc, "/label_code.R"))
    if(max(lesmask_dil_mask) > 0){
      label_result = label_lesion(lesmask_dil_mask, mimosa_reg_epi)
      writenii(label_result, paste0(reg.epi.out.dir,"/lesions_reg_epi_labeled"))
    }else{
      print(paste0("patient: ",p, " mimosa segmentation failed"))
    }
  }

  # WhiteStripe EPI Space
  white.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/whitestripe/EPI_space")
  #if(!file.exists(white.epi.out.dir)){
  # only create directory if it doesn't already exist - EAH 5/7/26
  if (!dir.exists(white.epi.out.dir)) {
    dir.create(white.epi.out.dir,showWarnings = FALSE)
  }
  
  # file check - EAH 5/6/26
  t1_ws_path = paste0(white.epi.out.dir,"/t1_n4_reg_epi_WS.nii.gz")
  flair_ws_path = paste0(white.epi.out.dir,"/flair_n4_reg_epi_WS.nii.gz")
  phase_n4_reg_brain_ws_path = paste0(white.epi.out.dir,"/phase_n4_WS_T2.nii.gz")
  
  files_ws_epi = c(t1_ws = t1_ws_path, flair_ws = flair_ws_path, phase_n4_reg_brain_ws = phase_n4_reg_brain_ws_path)
  missing_files_ws_epi = files_ws_epi[!file.exists(c(t1_ws_path, flair_ws_path, phase_n4_reg_brain_ws_path))]

  ## WhiteStripe T1, WS EPI using T1 indices
  ind = whitestripe(t1_reg_epi, "T1")
  
  if ('t1_ws' %in% names(missing_files_ws_epi)) {
    t1_ws = whitestripe_norm(t1_reg_epi, ind$whitestripe.ind)
    writenii(t1_ws,paste0(white.epi.out.dir,"/t1_n4_reg_epi_WS.nii.gz"))
  } else {
    t1_ws = readnii(t1_ws_path)
  }

  ## WhiteStripe FL, WS EPI using FL indices
  ind = whitestripe(flair_reg_epi, "T2")
  if ('flair_ws' %in% names(missing_files_ws_epi)) {
    flair_ws = whitestripe_norm(flair_reg_epi, ind$whitestripe.ind)
    writenii(flair_ws,paste0(white.epi.out.dir,"/flair_n4_reg_epi_WS.nii.gz"))
  } else {
    flair_ws = readnii(flair_ws_path)
  }

  ## WhiteStripe Phase using T2 indices
  ind2 = whitestripe(phase_n4_brain, "T2")
  if ('phase_n4_reg_brain_ws' %in% names(missing_files_ws_epi)) {
    phase_n4_reg_brain_ws = whitestripe_norm(phase_n4_brain, ind2$whitestripe.ind)
    writenii(phase_n4_reg_brain_ws, paste0(white.epi.out.dir, "/phase_n4_WS_T2.nii.gz"))
  } else {
    phase_n4_reg_brain_ws = readnii(phase_n4_reg_brain_ws_path)
  }
  #}

}else if(argv$step == "PRL_run"){
  my_path = paste0(argv$lesioncenter, "/")
  source_files = list.files(my_path)
  purrr::map(paste0(my_path, source_files), source)
  ## Find PRL
  p = argv$participant
  ses = argv$session
  prl.out.dir = paste0(main_path, "/data/", p, "/", ses, "/prl")
  # only create directory if it doesn't already exist - EAH 5/7/26
  if (!dir.exists(prl.out.dir)) {
    dir.create(prl.out.dir,showWarnings = FALSE)
  }

  # file check - EAH 5/7/26
  findprls_out_path = paste0(prl.out.dir,"/findprls_out_dil.rds")
  preds_path = paste0(prl.out.dir,"/",p,"_preds.csv")
  
  files_prl = c(findprls_out = findprls_out_path, preds = preds_path)
  missing_files_prl = files_prl[!file.exists(c(findprls_out_path, preds_path))]
  
  if ('findprls_out' %in% names(missing_files_prl)) {
    pretrainedmodel = readRDS(argv$aprlpath)
    source(paste0(argv$helpfunc, "/findprls_final.R")) 
    source(paste0(argv$helpfunc, "/extract_ria.R"))
  
    reg.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/EPI_space")
    white.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/whitestripe/EPI_space")
    label_result = readnii(paste0(reg.epi.out.dir,"/lesions_reg_epi_labeled"))
    findprls_out = findprls(lesmask = label_result, 
                            phasefile = paste0(white.epi.out.dir, "/phase_n4_WS_T2"),
                            pretrainedmodel = pretrainedmodel)
    saveRDS(findprls_out, paste0(prl.out.dir,"/findprls_out_dil.rds"))
  } else {
    findprls_out = readRDS(findprls_out_path)
  }
  
  if ('preds' %in% names(missing_files_prl)) {
    preds = findprls_out$preds
    write.csv(preds,paste0(prl.out.dir,"/",p,"_preds.csv"))
  } else {
    preds = read.csv(preds_path)
  }
}else if(argv$step == "consolidation"){
  if(!file.exists(paste0(main_path, "/stats"))){
    dir.create(paste0(main_path, "/stats"))
  }
  
  # only write file if it doesn't already exist - EAH 5/7/26
  if(!file.exists(paste0(main_path, "/stats/prl_probability.csv"))){
    prl_files = list.files(paste0(main_path, "/data"), pattern = "*_preds.csv", recursive = TRUE, full.names = TRUE) 
    prl_con = lapply(prl_files, function(x){
      sub_file = read_csv(x)
      subj = str_split(x, "/")[[1]][length(str_split(x, "/")[[1]]) - 2]
      sub_file = sub_file %>% mutate(subject = subj)
    }) %>% bind_rows()
    colnames(prl_con) = c("lesion_id", "rim_neg", "rim_pos", "subject")
    write_csv(prl_con, paste0(main_path, "/stats/prl_probability.csv"))
  }
}


