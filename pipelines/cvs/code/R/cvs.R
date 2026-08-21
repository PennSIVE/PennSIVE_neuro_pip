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


p <- arg_parser("Running Central Vein Sign (CVS) pipeline to obtain cvs probability.", hide.opts = FALSE)
p <- add_argument(p, "--mainpath", short = '-m', help = "Specify the main path where MRI images can be found.")
p <- add_argument(p, "--participant", short = '-p', help = "Specify the subject id.")
p <- add_argument(p, "--session", short = '-s', help = "Specify the session id.")
p <- add_argument(p, "--t1", help = "Specify the T1 sequence name.")
p <- add_argument(p, "--t2", help = "Specify the T2 sequence name.")
p <- add_argument(p, "--flair", help = "Specify the FLAIR sequence name.")
p <- add_argument(p, "--epi", help = "Specify the EPI sequence name.")
p <- add_argument(p, "--n4", help = "Specify whether to run bias correction step.", default = TRUE)
p <- add_argument(p, "--skullstripping", short = '-s', help = "Specify whether to run skull stripping step.", default = TRUE)
p <- add_argument(p, "--stype", help = "Specify which skullstripping method to use.", default = 'hdbet')
p <- add_argument(p, "--registration", short = '-r', help = "Specify whether to run registration step.", default = TRUE)
p <- add_argument(p, "--whitestripe", short = '-w', help = "Specify whether to run whitestripe step.", default = TRUE)
p <- add_argument(p, "--mimosa", help = "Specify whether to run mimosa segmentation step.", default = TRUE)
p <- add_argument(p, "--threshold", help = "Specify the threshold used to generate mimosa mask.", default = 0.2)
p <- add_argument(p, "--csf", help = "Specify whether to extract CSF mask.", default = TRUE)
p <- add_argument(p, "--step", help = "Specify the step of cvs pipeline. estimation or consolidation.", default = "estimation")
p <- add_argument(p, "--lesioncenter", help = "Provide the path to the lesioncenter package.")
p <- add_argument(p, "--mpath", help = "Specify the path to the trained mimosa model.")
p <- add_argument(p, "--helpfunc", help = "Specify the path to the help functions.")
argv <- parse_args(p)

# Read in Files
main_path = argv$mainpath
if (argv$step == "estimation"){
  # Load lesion center package
  my_path = paste0(argv$lesioncenter, "/")
  source_files = list.files(my_path)
  purrr::map(paste0(my_path, source_files), source)
  p = argv$participant
  ses = argv$session
  message('Checking inputs...')
  if(is.na(argv$t1)) stop("Missing T1 sequence!")else{
    t1 = readnii(paste0(main_path, "/data/", p,  "/", ses, "/anat/", argv$t1))
  }
  if(is.na(argv$flair)) stop("Missing FLAIR sequence!")else{
    flair = readnii(paste0(main_path, "/data/", p,  "/", ses, "/anat/", argv$flair))
  }

  if(is.na(argv$epi)) stop("Missing EPI sequence!")else{
    epi = read_rpi(paste0(main_path, "/data/", p,  "/", ses, "/anat/", argv$epi))
  }

  if(!is.na(argv$t2)){
    t2 = readnii(paste0(main_path, "/data/", p,  "/", ses, "/anat/", argv$t2))
  }

  # Bias Correction
  bias.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/bias_correction")
  if(argv$n4){
    message('Starting bias correction...')
    # only create directory if it doesn't already exist - EAH 5/6/26
    if (!dir.exists(bias.out.dir)) {
      dir.create(bias.out.dir,showWarnings = FALSE)
    }
    
    # file check - EAH 5/6/26
    t1_biascorrect_path = paste0(bias.out.dir,"/T1_n4.nii.gz")
    flair_biascorrect_path = paste0(bias.out.dir,"/FLAIR_n4.nii.gz")
    epi_biascorrect_path = paste0(bias.out.dir,"/EPI_n4.nii.gz")
    t2_biascorrect_path = paste0(bias.out.dir,"/T2_n4.nii.gz")
    
    if (is.na(argv$t2)) {
      files_n4 = c(t1_biascorrect = t1_biascorrect_path, flair_biascorrect = flair_biascorrect_path, epi_biascorrect = epi_biascorrect_path)
      missing_files_n4 = files_n4[!file.exists(c(t1_biascorrect_path, flair_biascorrect_path, epi_biascorrect_path))]
    } else {
      files_n4 = c(t1_biascorrect = t1_biascorrect_path, flair_biascorrect = flair_biascorrect_path, epi_biascorrect = epi_biascorrect_path, t2_biascorrect = t2_biascorrect_path)
      missing_files_n4 = files_n4[!file.exists(c(t1_biascorrect_path, flair_biascorrect_path, epi_biascorrect_path, t2_biascorrect_path))]
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
    
    if ('epi_biascorrect' %in% names(missing_files_n4)) {
      epi_biascorrect = bias_correct(file = epi,
                                   correction = "N4",
                                   verbose = TRUE)
      writenii(epi_biascorrect,paste0(bias.out.dir,"/EPI_n4.nii.gz"))
    } else {
      epi_biascorrect = readnii(epi_biascorrect_path)
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
    t1_biascorrect = readnii(paste0(main_path, "/data/", p,  "/", ses, "/bias_correction/T1_n4.nii.gz"))
    flair_biascorrect = readnii(paste0(main_path, "/data/", p,  "/", ses, "/bias_correction/FLAIR_n4.nii.gz"))
    if(!is.na(argv$t2)){
      t2_biascorrect = readnii(paste0(main_path, "/data/", p,  "/", ses, "/bias_correction/T2_n4.nii.gz")) 
    }
    existing_files = list.files(paste0(main_path, "/data/", p,  "/", ses, "/bias_correction"))
    epi_file = existing_files[which(grepl("EPI_n4*", existing_files))]
    if(length(epi_file) == 0){
      epi_biascorrect = bias_correct(file = epi,
                                 correction = "N4",
                                 verbose = TRUE)
    writenii(epi_biascorrect,paste0(bias.out.dir,"/EPI_n4.nii.gz"))
    }else{
      epi_biascorrect = read_rpi(paste0(main_path, "/data/", p,  "/", ses, "/bias_correction/EPI_n4.nii.gz"))
    }
  }

  message('Finished bias correction')


  # Skull Stripping
  brain.out.dir = paste0(main_path, "/data/", p, "/", ses, "/t1_brain")

  if (!argv$skullstripping) {
    # read pre-existing skull-stripped brain
    brain_paths = list.files(brain.out.dir, recursive = TRUE, full.names = TRUE)

    # Match *brain.nii.gz but exclude mask files
    brain_path = brain_paths[grepl("brain\\.nii\\.gz$", brain_paths) &
                              !grepl("(brainmask|brain_mask)\\.nii\\.gz$", brain_paths)]

    # Accept both FSL BET convention (T1_brainmask.nii.gz)
    #        and HD-BET convention  (T1_brain_mask.nii.gz)
    brain_mask_path = brain_paths[grepl("(brainmask|brain_mask)\\.nii\\.gz$", brain_paths)]

    if (length(brain_path) == 0) {
      stop("No brain file found in: ", brain.out.dir,
          "\nExpected a file matching *brain.nii.gz (e.g. T1_brain.nii.gz)")
    }
    if (length(brain_mask_path) == 0) {
      stop("No brain mask found in: ", brain.out.dir,
          "\nExpected T1_brainmask.nii.gz (FSL BET) or T1_brain_mask.nii.gz (HD-BET)")
    }

    t1_brain = readnii(brain_path)
    brain_mask       = readnii(brain_mask_path)
    t1_brain = bias_correct(file = t1_brain, correction = "N4", verbose = TRUE)

    if (!file.exists(paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))) {
      writenii(t1_brain, paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))
    }
  }

  if (argv$skullstripping) {
    message('Starting skullstripping...')
    # run skull stripping
    t1_brain_path   = paste0(brain.out.dir, "/T1_brain.nii.gz")
    brain_mask_path = paste0(brain.out.dir, "/T1_brainmask.nii.gz")

    files_skullstripping = c(t1_brain   = t1_brain_path,
                              brain_mask = brain_mask_path)
    missing_files_skullstripping = files_skullstripping[
      !file.exists(c(t1_brain_path, brain_mask_path))
    ]

    if (!dir.exists(brain.out.dir)) {
      dir.create(brain.out.dir, showWarnings = FALSE)
    }

    # FSL BET
    if (argv$stype == "fslbet") {
      message("Running FSL BET skull stripping on bias-corrected T1...")

      if ('t1_brain' %in% names(missing_files_skullstripping)) {
        t1_brain = fslbet_robust(t1_biascorrect, reorient = FALSE, correct = FALSE)
        writenii(t1_brain, t1_brain_path)
        writenii(t1_brain, paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))
      } else {
        t1_brain = readnii(t1_brain_path)
      }

      if ('brain_mask' %in% names(missing_files_skullstripping)) {
        brain_mask = t1_brain > 0
        writenii(brain_mask, brain_mask_path)
      } else {
        brain_mask = readnii(brain_mask_path)
      }

    # HD-BET
    } else if (argv$stype == "hdbet") {

      needs_hdbet = ('t1_brain'   %in% names(missing_files_skullstripping)) ||
                    ('brain_mask' %in% names(missing_files_skullstripping))

      if (needs_hdbet) {
        t1_input = paste0(bias.out.dir, "/T1_n4.nii.gz")

        if (!file.exists(t1_input)) {
          stop("Bias-corrected T1 not found at: ", t1_input,
              "\nEnsure --n4 TRUE is set so bias correction runs before skull stripping.")
        }

        message("Running HD-BET skull stripping on bias-corrected T1...")
        ret = system(paste("hd-bet",
                    "-i", t1_input,
                    "-o", t1_brain_path,
                    "-device cpu --disable_tta --save_bet_mask"))
        if (ret != 0) stop("HD-BET exited with non-zero status: ", ret)

        # HD-BET writes the mask as <output>_bet.nii.gz; rename to match PennSIVE_neuro_pip brainmask naming convention
        hdbet_mask_raw = sub("\\.nii\\.gz$", "_bet.nii.gz", t1_brain_path)
        if (file.exists(hdbet_mask_raw) && !file.exists(brain_mask_path)) {
          file.rename(hdbet_mask_raw, brain_mask_path)
        }

        if (!file.exists(paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))) {
          writenii(readnii(t1_brain_path), paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))
        }
      }

      t1_brain = readnii(t1_brain_path)
      brain_mask       = readnii(brain_mask_path)

    } else {
      stop("Unknown --stype '", argv$stype, "'. Use 'fslbet' or 'hdbet'.")
    }
  }

  message('Finished skullstripping')



  # Registration to FLAIR Space
  reg.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/registration/FLAIR_space")
  if (argv$registration){
    message('Starting registration...')
    # only create directory if it doesn't already exist - EAH 5/6/26
    if (!dir.exists(reg.out.dir)) {
      dir.create(reg.out.dir,showWarnings = FALSE, recursive = TRUE)
    }
    
    # file check - EAH 5/6/26
    t1_reg_path = paste0(reg.out.dir,"/t1_n4_brain_reg_flair.nii.gz")
    brainmask_reg_path = paste0(reg.out.dir,"/brainmask_reg_flair.nii.gz")
    flair_n4_brain_path = paste0(reg.out.dir,"/flair_n4_brain.nii.gz")
    t2_n4_brain_path = paste0(reg.out.dir,"/t2_n4_brain_reg_flair.nii.gz")
    
    if (is.na(argv$t2)) {
      files_reg_flair = c(t1_reg = t1_reg_path, brainmask_reg = brainmask_reg_path, flair_n4_brain = flair_n4_brain_path)
      missing_files_reg_flair = files_reg_flair[!file.exists(c(t1_reg_path, brainmask_reg_path, flair_n4_brain_path))]
    } else {
      files_reg_flair = c(t1_reg = t1_reg_path, brainmask_reg = brainmask_reg_path, flair_n4_brain = flair_n4_brain_path, t2_n4_brain = t2_n4_brain_path)
      missing_files_reg_flair = files_reg_flair[!file.exists(c(t1_reg_path, brainmask_reg_path, flair_n4_brain_path, t2_n4_brain_path))]
    }
    

    ## Register T1 to FLAIR space 
    
    t1_to_flair = registration(filename = t1_biascorrect,
                             template.file = flair_biascorrect,
                             typeofTransform = "Rigid", remove.warp = FALSE,
                             outprefix=paste0(reg.out.dir,"/t1_reg_to_flair")) 
    if ('t1_reg' %in% names(missing_files_reg_flair)) {
      t1_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(t1_brain),
                                          transformlist = t1_to_flair$fwdtransforms, interpolator = "welchWindowedSinc"))
      writenii(t1_reg, paste0(reg.out.dir,"/t1_n4_brain_reg_flair"))
    } else {
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
      if ('t2_n4_brain' %in% names(missing_files_reg_flair)) {
        t2_to_flair = registration(filename = t2_biascorrect,
                                 template.file = flair_biascorrect,
                                 typeofTransform = "Rigid", remove.warp = FALSE,
                                 outprefix=paste0(reg.out.dir,"/t2_reg_to_flair"))
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

  message('Finished registration')
  

  # WhiteStripe normalize data
  white.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/whitestripe/FLAIR_space")
  if(argv$whitestripe){
    message('Starting whitestripe...')
    # only create directory if it doesn't already exist - EAH 5/6/26
    if (!dir.exists(white.out.dir)) {
      dir.create(white.out.dir,showWarnings = FALSE, recursive = TRUE)
    }
    
    # file check - EAH 5/6/26
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

    message('Finished whitestripe')

  # Mimosa
  mim.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/mimosa")
  if(argv$mimosa){
    message('Starting MIMoSA...')
    # only create directory if it doesn't already exist - EAH 5/6/26
    if (!dir.exists(mim.out.dir)) {
      dir.create(mim.out.dir,showWarnings = FALSE)
    }

    # file check - EAH 5/6/26
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

  message('Finished MIMoSA')

  # Extract CSF
  csf.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/CSF")
  if(argv$csf){
    # only create directory if it doesn't already exist - EAH 5/6/26
    if (!dir.exists(csf.out.dir)) {
      dir.create(csf.out.dir,showWarnings = FALSE)
    }
    
    # file check - EAH 5/6/26
    csf_path = paste0(csf.out.dir,"/csf.nii.gz")

    files_csf = c(csf = csf_path)
    missing_files_csf = files_csf[!file.exists(c(csf_path))]
    
    if ('csf' %in% names(missing_files_csf)) {
      csf=fast(t1_reg,opts='--nobias')
      csf[csf!=1]= 0
      csf=ants2oro(labelClusters(oro2ants(csf),minClusterSize=300))
      csf[csf>0]= 1
      csf=(csf!=1)
      csf=fslerode(csf, kopts = paste("-kernel boxv",2), verbose = TRUE)
      csf=(csf==0)
      writenii(csf, paste0(csf.out.dir,"/csf"))
    } else {
      csf = readnii(csf_path)
    }
  }else{
    csf = readnii(paste0(csf.out.dir,"/csf"))
  }

  # Split Confluent Lesions
  les = lesioncenters(probmap, probmap>argv$threshold, parallel=FALSE, cores=1, c3d=F)$lesioncenters
  lables=ants2oro(labelClusters(oro2ants(les>0),minClusterSize=27))
  for(j in 1:max(lables)){
    if(sum(csf[lables==j])>0){
      lables[lables==j] = 0
      }
    }
  les = lables > 0

  # Register to EPI Space
  reg.epi.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/registration/EPI_space")
  
  # only create directory if it doesn't already exist - EAH 5/6/26
  if (!dir.exists(reg.epi.out.dir)) {
    dir.create(reg.epi.out.dir,showWarnings = FALSE)
  }
  
  # file check - EAH 5/6/26
  
  brainmask_reg_epi_path = paste0(reg.epi.out.dir,'/brainmask_reg_epi.nii.gz')
  epi_n4_brain_path = paste0(reg.epi.out.dir,'/epi_n4_brain.nii.gz')
  t1_reg_epi_path = paste0(reg.epi.out.dir, "/t1_reg_epi.nii.gz")
  flair_reg_epi_path = paste0(reg.epi.out.dir, "/flair_reg_epi.nii.gz")
  mimosa_reg_epi_path = paste0(reg.epi.out.dir, "/mimosa_reg_epi.nii.gz")
  mimosa_mask_reg_epi_path = paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi.nii.gz")
  les_reg_epi_path = paste0(reg.epi.out.dir, "/les_reg_epi.nii.gz")
  
  files_reg_epi = c(brainmask_reg_epi = brainmask_reg_epi_path, epi_n4_brain = epi_n4_brain_path, t1_reg_epi = t1_reg_epi_path, flair_reg_epi = flair_reg_epi_path,
                    mimosa_reg_epi = mimosa_reg_epi_path, mimosa_mask_reg_epi = mimosa_mask_reg_epi_path, les_reg_epi = les_reg_epi_path)
  missing_files_reg_epi = files_reg_epi[!file.exists(c(brainmask_reg_epi_path, epi_n4_brain_path, t1_reg_epi_path, flair_reg_epi_path,
                                                       mimosa_reg_epi_path, mimosa_mask_reg_epi_path, les_reg_epi_path))]

  
  flair_to_epi = registration(filename = flair_biascorrect,
                                template.file = epi_biascorrect,
                                typeofTransform = "Rigid", remove.warp = FALSE) ### rigid

  if ('brainmask_reg_epi' %in% names(missing_files_reg_epi)) {
    brainmask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_biascorrect), moving = oro2ants(brainmask_reg),
                                                  transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
    writenii(brainmask_reg_epi, paste0(reg.epi.out.dir,'/brainmask_reg_epi'))
  } else {
    brainmask_reg_epi = readnii(brainmask_reg_epi_path)
  }
  
  if ('epi_n4_brain' %in% names(missing_files_reg_epi)) {
    epi_n4_brain = epi_biascorrect * brainmask_reg_epi
    writenii(epi_n4_brain, paste0(reg.epi.out.dir,'/epi_n4_brain'))
  } else {
    epi_n4_brain = readnii(epi_n4_brain_path)
  }

  if ('t1_reg_epi' %in% names(missing_files_reg_epi)) {
    t1_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(t1_reg),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
    writenii(t1_reg_epi, paste0(reg.epi.out.dir, "/t1_reg_epi"))
  } else {
    t1_reg_epi = readnii(t1_reg_epi_path)
  }
  
  if ('flair_reg_epi' %in% names(missing_files_reg_epi)) {
    flair_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(flair_n4_brain),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
    writenii(flair_reg_epi, paste0(reg.epi.out.dir, "/flair_reg_epi"))
  } else {
    flair_reg_epi = readnii(flair_reg_epi_path)
  }
  
  if ('mimosa_reg_epi' %in% names(missing_files_reg_epi)) {
    mimosa_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(probmap),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
    writenii(mimosa_reg_epi, paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
  } else {
    mimosa_reg_epi = readnii(mimosa_reg_epi_path)
  }
  
  if ('mimosa_mask_reg_epi' %in% names(missing_files_reg_epi)) {
    mimosa_mask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(probmap > argv$threshold),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
    writenii(mimosa_mask_reg_epi, paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
  } else {
    mimosa_mask_reg_epi = readnii(mimosa_mask_reg_epi_path)
  }
  
  if ('les_reg_epi' %in% names(missing_files_reg_epi)) {
    les_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(epi_n4_brain), moving = oro2ants(les),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
    writenii(les_reg_epi, paste0(reg.epi.out.dir, "/les_reg_epi"))
  } else {
    les_reg_epi = readnii(les_reg_epi_path)
  }

  # CVS Score Calculation
  message('Starting CVS calculation...')
  cvs.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/cvs")
  # only create directory if it doesn't already exist - EAH 5/6/26
  if (!dir.exists(cvs.out.dir)) {
    dir.create(cvs.out.dir,showWarnings = FALSE)
  }
  
  # file check - EAH 5/6/26
  frangi_path = paste0(cvs.out.dir,"/frangi.nii.gz")
  dtb_path = paste0(cvs.out.dir,"/dtb.nii.gz")
  lables_path = paste0(cvs.out.dir,"/lables.nii.gz")
  probles_path = paste0(cvs.out.dir,"/cvs_probmap.nii.gz")
  summary_df_path = paste0(cvs.out.dir,"/cvs_biomarker.csv")
  summary_df_lesion_path = paste0(cvs.out.dir,"/cvs_biomarker_lesion.csv")
  
  files_cvs = c(frangi = frangi_path, dtb = dtb_path, lables = lables_path, probles = probles_path, summary_df = summary_df_path, summary_df_lesion = summary_df_lesion_path)
  missing_files_cvs = files_cvs[!file.exists(c(frangi_path, dtb_path, lables_path, probles_path, summary_df_path, summary_df_lesion_path))]
  
  ## Obtain vein map 
  source(paste0(argv$helpfunc, "/helperfunctions.R"))
  if ('frangi' %in% names(missing_files_cvs)) {
    frangi = frangifilternoc3d(image = epi_biascorrect, mask = brainmask_reg_epi)
    frangi[frangi<0]= 0
    writenii(frangi, paste0(cvs.out.dir,"/frangi"))
  } else {
    frangi = readnii(frangi_path)
  }

  ## Obtain distance-to-the-boundary map 
  if ('dtb' %in% names(missing_files_cvs)) {
    dtb=dtboundary(les_reg_epi)
    writenii(dtb, paste0(cvs.out.dir,"/dtb"))
  } else {
    dtb = readnii(dtb_path)
  }
  
  ## Perform permutation procedure to get CVS probabilities
  if ('lables' %in% names(missing_files_cvs)) {
    lables = ants2oro(labelClusters(oro2ants(les_reg_epi),minClusterSize=27))
    writenii(lables, paste0(cvs.out.dir,"/lables"))
  } else {
    lables = readnii(lables_path)
  }
  
  if ('probles' %in% names(missing_files_cvs)) {
    probles = lables
    avprob = NULL
    maxles = max(as.vector(lables))
    for(j in 1:maxles){
      # get true coherence for lesion j
      frangsub = frangi[lables==j]
      centsub = dtb[lables==j]
      coords = which(lables==j, arr.ind=T)
      prod = frangsub*centsub
      score = sum(prod)
      # get 1000 null coherence values for lesion j
      nullscores = as.vector(unlist(lapply(1:1000, getnulldist, centsub, coords, frangsub)))
      # get CVS probability for lesion j
      lesprob = sum(nullscores<score)/length(nullscores)
      avprob = c(avprob,lesprob)
      probles[lables==j]= lesprob
      print(paste0("Done with lesion ", j ," of ", maxles))
    }
    writenii(probles, paste0(cvs.out.dir,"/cvs_probmap.nii.gz"))
  } else {
    probles = readnii(probles_path)
  }
  
  #saveRDS(avprob, paste0(cvs.out.dir,"/cvs_avprob"))
  
  if ('summary_df' %in% names(missing_files_cvs)) {
    cvs.biomarker = mean(avprob)
    cvs.biomarker.lesion = avprob
    numles = maxles
    lesion_id = 1:numles
    subject_id = p
    session_id = ses
    summary_df = data.frame(cbind(subject_id, session_id, numles, cvs.biomarker))
    write_csv(summary_df, paste0(cvs.out.dir,"/cvs_biomarker.csv"))
  } else {
    summary_df = read.csv(summary_df_path)
  }
  
  if ('summary_df_lesion' %in% names(missing_files_cvs)) {
    summary_df_lesion = data.frame(cbind(subject_id = rep(subject_id, numles), session_id = rep(session_id, numles), lesion_id, cvs.score = cvs.biomarker.lesion))
    write_csv(summary_df_lesion, paste0(cvs.out.dir,"/cvs_biomarker_lesion.csv"))
  } else {
    summary_df_lesion = read.csv(summary_df_lesion_path)
  }

  message('Finished CVS calculation')

  # WhiteStripe EPI Space
  white.epi.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/whitestripe/EPI_space")
  
  # only create directory if it doesn't already exist - EAH 5/6/26
  if (!dir.exists(white.epi.out.dir)) {
    dir.create(white.epi.out.dir,showWarnings = FALSE)
  }
  
  # file check - EAH 5/6/26
  t1_ws_path = paste0(white.epi.out.dir,"/t1_n4_reg_epi_WS.nii.gz")
  epi_ws_T1_path = paste0(white.epi.out.dir,"/epi_n4_WS_T1.nii.gz")
  flair_ws_path = paste0(white.epi.out.dir,"/flair_n4_reg_epi_WS.nii.gz")
  epi_ws_T2_path = paste0(white.epi.out.dir,"/epi_n4_WS_FL.nii.gz")
  
  files_ws_epi = c(t1_ws = t1_ws_path, epi_ws_T1 = epi_ws_T1_path, flair_ws = flair_ws_path, epi_ws_T2 = epi_ws_T2_path)
  missing_files_ws_epi = files_ws_epi[!file.exists(c(t1_ws_path, epi_ws_T1_path, flair_ws_path, epi_ws_T2_path))]

  ## WhiteStripe T1, WS EPI using T1 indices
  ind = whitestripe(t1_reg_epi, "T1")
  
  if ('t1_ws' %in% names(missing_files_ws_epi)) {
    t1_ws = whitestripe_norm(t1_reg_epi, ind$whitestripe.ind)
    writenii(t1_ws,paste0(white.epi.out.dir,"/t1_n4_reg_epi_WS.nii.gz"))
  } else {
    t1_ws = readnii(t1_ws_path)
  }
  
  if ('epi_ws_T1' %in% names(missing_files_ws_epi)) {
    epi_ws_T1 = whitestripe_norm(epi_n4_brain, ind$whitestripe.ind)
    writenii(epi_ws_T1,paste0(white.epi.out.dir,"/epi_n4_WS_T1.nii.gz"))
  } else {
    epi_ws_T1 = readnii(epi_ws_T1_path)
  }
  
  ## WhiteStripe FL, WS EPI using FL indices
  ind = whitestripe(flair_reg_epi, "T2")
  
  if ('flair_ws' %in% names(missing_files_ws_epi)) {
    flair_ws = whitestripe_norm(flair_reg_epi, ind$whitestripe.ind)
    writenii(flair_ws,paste0(white.epi.out.dir,"/flair_n4_reg_epi_WS.nii.gz"))
  } else {
    flair_ws = readnii(flair_ws_path)
  }
  
  if ('epi_ws_T2' %in% names(missing_files_ws_epi)) {
    epi_ws_T2 = whitestripe_norm(epi_n4_brain, ind$whitestripe.ind)
    writenii(epi_ws_T2,paste0(white.epi.out.dir,"/epi_n4_WS_FL.nii.gz"))
  } else {
    epi_ws_T2 = readnii(epi_ws_T2_path)
  }
}else if(argv$step == "consolidation"){
  if(!file.exists(paste0(main_path, "/stats"))){
    dir.create(paste0(main_path, "/stats"))
  }
  
  # only write file if it doesn't already exist - EAH 5/6/26
  if(!file.exists(paste0(main_path, "/stats/cvs_score.csv"))){
    cvs_con = list.files(paste0(main_path, "/data"), pattern = "cvs_biomarker.csv", recursive = TRUE, full.names = TRUE) %>% read_csv() %>% bind_rows()
    write_csv(cvs_con, paste0(main_path, "/stats/cvs_score.csv"))
  }
}