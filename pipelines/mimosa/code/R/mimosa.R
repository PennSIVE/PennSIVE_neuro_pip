.libPaths(c("/misc/appl/R-4.5/lib64/R/library", .libPaths()))

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

p <- arg_parser("Running MIMoSA Model to segement White Matter MS lesions", hide.opts = FALSE)
p <- add_argument(p, "--mainpath", short = '-m', help = "Specify the main path where MRI images can be found.")
p <- add_argument(p, "--participant", short = '-p', help = "Specify the subject id.")
p <- add_argument(p, "--session", short = '-s', help = "Specify the session.")
p <- add_argument(p, "--t1", help = "Specify the T1 sequence name")
p <- add_argument(p, "--t2", help = "Specify the T2 sequence name")
p <- add_argument(p, "--flair", help = "Specify the FLAIR sequence name")
p <- add_argument(p, "--n4", help = "Specify whether to run bias correction step.", default = TRUE)
p <- add_argument(p, "--skullstripping", help = "Specify whether to run skull stripping step.", default = TRUE)
p <- add_argument(p, "--stype", help = "Specify which skullstripping method to use.", default = 'hdbet')
p <- add_argument(p, "--registration", short = '-r', help = "Specify whether to run registration step.", default = TRUE)
p <- add_argument(p, "--whitestripe", short = '-w', help = "Specify whether to run whitestripe step.", default = TRUE)
p <- add_argument(p, "--threshold", help = "Specify the threshold used to generate mimosa mask.", default = 0.2)
p <- add_argument(p, "--container", help = "Specify the container to use: singularity, docker, local, cluster.", default = 'cluster')
p <- add_argument(p, "--mpath", help = "Specify the path to the trained mimosa model.")
p <- add_argument(p, "--toolpath", help = "Specify the path to the saved pipeline folder, eg: /path/to/folder")
argv <- parse_args(p)

# Read in Files
main_path = argv$mainpath
p = argv$participant
s = argv$session
model_path = argv$mpath
tool_path = argv$toolpath
message('Checking inputs...')
if(is.na(argv$t1)) stop("Missing T1 sequence!")else{
  t1 = readnii(paste0(main_path, "/data/", p, "/", s, "/anat/", argv$t1))
}
if(is.na(argv$flair)) stop("Missing FLAIR sequence!")else{
  flair = readnii(paste0(main_path, "/data/", p, "/", s, "/anat/", argv$flair))
}

if(!is.na(argv$t2)){
  t2 = readnii(paste0(main_path, "/data/", p, "/", s, "/anat/", argv$t2))
}

# Bias Correction
if(argv$n4){
  message('Starting bias correction...')
  bias.out.dir = paste0(main_path, "/data/", p, "/", s, "/bias_correction")
  # only create directory if it doesn't already exist
  if (!dir.exists(bias.out.dir)) {
    dir.create(bias.out.dir,showWarnings = FALSE)
  }
  
  # file check
  t1_biascorrect_path = paste0(bias.out.dir,"/T1_n4.nii.gz")
  flair_biascorrect_path = paste0(bias.out.dir,"/FLAIR_n4.nii.gz")
  t2_biascorrect_path = paste0(bias.out.dir,"/T2_n4.nii.gz")
  
  if (is.na(argv$t2)) {
    files_n4 = c(t1_biascorrect = t1_biascorrect_path, flair_biascorrect = flair_biascorrect_path)
    missing_files_n4 = files_n4[!file.exists(c(t1_biascorrect_path, flair_biascorrect_path))]
  } else {
    files_n4 = c(t1_biascorrect = t1_biascorrect_path, flair_biascorrect = flair_biascorrect_path, t2_biascorrect = t2_biascorrect_path)
    missing_files_n4 = files_n4[!file.exists(c(t1_biascorrect_path, flair_biascorrect_path, t2_biascorrect_path))]
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
  bias.out.dir = paste0(main_path, "/data/", p, "/", s, "/bias_correction")
  t1_biascorrect = readnii(paste0(main_path, "/data/", p, "/", s, "/bias_correction/T1_n4.nii.gz"))
  flair_biascorrect = readnii(paste0(main_path, "/data/", p, "/", s, "/bias_correction/FLAIR_n4.nii.gz"))
  if(!is.na(argv$t2)){
    t2_biascorrect = readnii(paste0(main_path, "/data/", p, "/", s, "/bias_correction/T2_n4.nii.gz")) 
  }
}

message('Finished bias correction')


# Skull Stripping
brain.out.dir = paste0(main_path, "/data/", p, "/", s, "/t1_brain")

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
reg.out.dir = paste0(main_path, "/data/", p, "/", s, "/registration/FLAIR_space")

if (argv$registration) {
  message('Starting registration...')
  if (!dir.exists(reg.out.dir)) {
    dir.create(reg.out.dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Define output file paths
  t1_to_flair_mat    = paste0(reg.out.dir, "/t1_reg_to_flair0GenericAffine.mat")
  t1_reg_path        = paste0(reg.out.dir, "/t1_n4_brain_reg_flair.nii.gz")
  brainmask_reg_path = paste0(reg.out.dir, "/brainmask_reg_flair.nii.gz")
  flair_n4_brain_path = paste0(reg.out.dir, "/flair_n4_brain.nii.gz")
  t2_to_flair_mat    = paste0(reg.out.dir, "/t2_reg_to_flair0GenericAffine.mat")
  t2_n4_brain_path   = paste0(reg.out.dir, "/t2_n4_brain_reg_flair.nii.gz")

  # ── Register T1 to FLAIR space ──────────────────────────────────────────────
  # Run registration only if the .mat file doesn't exist yet.
  # We don't store the return object — the .mat file path is used directly
  # for antsApplyTransforms, which is more robust across ANTsR/extrantsr versions.
  if (!file.exists(t1_to_flair_mat)) {
    antsRegistration(
      fixed          = oro2ants(flair_biascorrect),
      moving         = oro2ants(t1_biascorrect),
      typeofTransform = "Rigid",
      outprefix      = paste0(reg.out.dir, "/t1_reg_to_flair")
    )
  }

  if (!file.exists(t1_reg_path)) {
    t1_reg = ants2oro(antsApplyTransforms(
      fixed         = oro2ants(flair_biascorrect),
      moving        = oro2ants(t1_brain),
      transformlist = t1_to_flair_mat,
      interpolator  = "welchWindowedSinc"
    ))
    writenii(t1_reg, paste0(reg.out.dir, "/t1_n4_brain_reg_flair"))
  } else {
    t1_reg = readnii(t1_reg_path)
  }

  if (!file.exists(brainmask_reg_path)) {
    brainmask_reg = ants2oro(antsApplyTransforms(
      fixed         = oro2ants(flair_biascorrect),
      moving        = oro2ants(brain_mask),
      transformlist = t1_to_flair_mat,
      interpolator  = "nearestNeighbor"
    ))
    writenii(brainmask_reg, paste0(reg.out.dir, "/brainmask_reg_flair"))
  } else {
    brainmask_reg = readnii(brainmask_reg_path)
  }

  if (!file.exists(flair_n4_brain_path)) {
    flair_n4_brain = flair_biascorrect
    flair_n4_brain[brainmask_reg == 0] = 0
    writenii(flair_n4_brain, paste0(reg.out.dir, "/flair_n4_brain"))
  } else {
    flair_n4_brain = readnii(flair_n4_brain_path)
  }

  # ── Register T2 to FLAIR space (optional) ───────────────────────────────────
  if (!is.na(argv$t2)) {

    if (!file.exists(t2_to_flair_mat)) {
      antsRegistration(
        fixed           = oro2ants(flair_biascorrect),
        moving          = oro2ants(t2_biascorrect),
        typeofTransform = "Rigid",
        outprefix       = paste0(reg.out.dir, "/t2_reg_to_flair")
      )
    }

    if (!file.exists(t2_n4_brain_path)) {
      t2_n4_brain = ants2oro(antsApplyTransforms(
        fixed         = oro2ants(flair_biascorrect),
        moving        = oro2ants(t2_biascorrect),
        transformlist = t2_to_flair_mat,
        interpolator  = "welchWindowedSinc"
      ))
      t2_n4_brain[brainmask_reg == 0] = 0
      writenii(t2_n4_brain, paste0(reg.out.dir, "/t2_n4_brain_reg_flair"))
    } else {
      t2_n4_brain = readnii(t2_n4_brain_path)
    }

  }

} else {
  t1_reg        = readnii(paste0(reg.out.dir, "/t1_n4_brain_reg_flair.nii.gz"))
  flair_n4_brain = readnii(paste0(reg.out.dir, "/flair_n4_brain.nii.gz"))
  brainmask_reg  = readnii(paste0(reg.out.dir, "/brainmask_reg_flair.nii.gz"))
  if (!is.na(argv$t2)) {
    t2_n4_brain = readnii(paste0(reg.out.dir, "/t2_n4_brain_reg_flair.nii.gz"))
  }
}

message('Finished registration')


# WhiteStripe normalize data
white.out.dir = paste0(main_path, "/data/", p, "/", s, "/whitestripe/FLAIR_space")
if(argv$whitestripe){
  message('Starting whitestripe...')
  # only create directory if it doesn't already exist
  if (!dir.exists(white.out.dir)) {
    dir.create(white.out.dir,showWarnings = FALSE, recursive = TRUE)
  }
  
  # file check
  t1_n4_brain_reg_ws_path = paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws.nii.gz")
  t2_n4_brain_reg_ws_path = paste0(white.out.dir, "/t2_n4_brain_reg_flair_ws.nii.gz")
  flair_n4_brain_ws_path = paste0(white.out.dir, "/flair_n4_brain_ws.nii.gz")
  
  if (is.na(argv$t2)) {
    files_ws_flair = c(t1_n4_brain_reg_ws = t1_n4_brain_reg_ws_path, flair_n4_brain_ws = flair_n4_brain_ws_path)
    missing_files_ws_flair = files_ws_flair[!file.exists(c(t1_n4_brain_reg_ws_path, flair_n4_brain_ws_path))]
  } else {
    files_ws_flair = c(t1_n4_brain_reg_ws = t1_n4_brain_reg_ws_path, t2_n4_brain_reg_ws = t2_n4_brain_reg_ws_path, flair_n4_brain_ws = flair_n4_brain_ws_path)
    missing_files_ws_flair = files_ws_flair[!file.exists(c(t1_n4_brain_reg_ws_path, t2_n4_brain_reg_ws_path, flair_n4_brain_ws_path))]
  }
  
  
  if ('t1_n4_brain_reg_ws' %in% names(missing_files_ws_flair)) {
    ind1 = whitestripe(t1_reg, "T1")
    t1_n4_brain_reg_ws = whitestripe_norm(t1_reg, ind1$whitestripe.ind)
    writenii(t1_n4_brain_reg_ws, paste0(white.out.dir,"/t1_n4_brain_reg_flair_ws"))
  } else {
    t1_n4_brain_reg_ws = readnii(t1_n4_brain_reg_ws_path)
  }
  
  if(!is.na(argv$t2)){
    if ('t2_n4_brain_reg_ws' %in% names(missing_files_ws_flair)) {
      ind2 = whitestripe(t2_n4_brain, "T2")
      t2_n4_brain_reg_ws = whitestripe_norm(t2_n4_brain, ind2$whitestripe.ind)
      writenii(t2_n4_brain_reg_ws, paste0(white.out.dir,"/t2_n4_brain_reg_flair_ws"))
    } else {
      t2_n4_brain_reg_ws = readnii(t2_n4_brain_reg_ws_path)
    }
  }
  
  if ('flair_n4_brain_ws' %in% names(missing_files_ws_flair)) {
    ind3 = whitestripe(flair_n4_brain, "T2")
    flair_n4_brain_ws = whitestripe_norm(flair_n4_brain, ind3$whitestripe.ind)
    writenii(flair_n4_brain_ws, paste0(white.out.dir,"/flair_n4_brain_ws"))
  } else {
    flair_n4_brain_ws = readnii(flair_n4_brain_ws_path)
  }
} else{
  t1_n4_brain_reg_ws = readnii(paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws"))
  if(!is.na(argv$t2)){
    t2_n4_brain_reg_ws = readnii(paste0(white.out.dir, "/t2_n4_brain_reg_flair_ws"))
  }
  flair_n4_brain_ws = readnii(paste0(white.out.dir, "/flair_n4_brain_ws"))
}

message('Finished whitestripe')


# Mimosa
mim.out.dir = paste0(main_path, "/data/", p, "/", s, "/mimosa")
message('Starting MIMoSA...')

# only create directory if it doesn't already exist
if (!dir.exists(mim.out.dir)) {
  dir.create(mim.out.dir,showWarnings = FALSE)
}

# file check
probmap_path = paste0(mim.out.dir,"/mimosa_prob")
mimosa_mask_path = paste0(mim.out.dir,"/mimosa_mask")

files_mimosa = c(probmap = probmap_path, mimosa_mask = mimosa_mask_path)
missing_files_mimosa = files_mimosa[!file.exists(c(probmap_path, mimosa_mask_path))]

if ('probmap' %in% names(missing_files_mimosa)) {
  mimosa = mimosa_data(brain_mask=brainmask_reg, FLAIR=flair_n4_brain_ws, T1=t1_n4_brain_reg_ws, gold_standard=NULL, normalize="no", cores = 1, verbose = TRUE)
  mimosa_df = mimosa$mimosa_dataframe
  cand_voxels = mimosa$top_voxels
  tissue_mask = mimosa$tissue_mask
  load(model_path) 
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

message('Finished MIMoSA')
