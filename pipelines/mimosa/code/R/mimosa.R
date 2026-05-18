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
p <- add_argument(p, "--skullstripping", help = "Specify whether to run skull stripping step.", default = FALSE)
p <- add_argument(p, "--registration", short = '-r', help = "Specify whether to run registration step.", default = TRUE)
p <- add_argument(p, "--whitestripe", short = '-w', help = "Specify whether to run whitestripe step.", default = TRUE)
p <- add_argument(p, "--threshold", help = "Specify the threshold used to generate mimosa mask.", default = 0.2)
p <- add_argument(p, "--mpath", help = "Specify the path to the trained mimosa model.")
argv <- parse_args(p)

# Read in Files
main_path = argv$mainpath
p = argv$participant
s = argv$session
model_path = argv$mpath
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
  
  bias.out.dir = paste0(main_path, "/data/", p, "/", s, "/bias_correction")
  # only create directory if it doesn't already exist - EAH 5/5/26
  if (!dir.exists(bias.out.dir)) {
    dir.create(bias.out.dir,showWarnings = FALSE)
  }
  
  # file check - EAH 5/5/26
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

# Skull Stripping
brain.out.dir = paste0(main_path, "/data/", p, "/", s, "/t1_brain")
if(!argv$skullstripping){
  brain_paths = list.files(paste0(main_path, "/data/", p, "/", s, "/t1_brain"), recursive = TRUE, full.names = TRUE)
  brain_path = brain_paths[which(grepl("*brain.nii.gz$", brain_paths))]
  brain_mask_path = brain_paths[which(grepl("*brainmask.nii.gz$", brain_paths))]
  t1_fslbet_robust = readnii(brain_path)
  brain_mask = readnii(brain_mask_path) 
  t1_fslbet_robust = bias_correct(file = t1_fslbet_robust,
                             correction = "N4",
                             verbose = TRUE)
  # only write file if it doesn't exist - EAH 5/6/26
  if (!file.exists(paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))) {
    writenii(t1_fslbet_robust,paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
  }
}

if (argv$skullstripping){
  # file check - EAH 5/5/26
  t1_fslbet_robust_path = paste0(brain.out.dir,"/T1_brain.nii.gz")
  brain_mask_path = paste0(brain.out.dir,"/T1_brainmask.nii.gz")
  
  files_skullstripping = c(t1_fslbet_robust = t1_fslbet_robust_path, brain_mask = brain_mask_path)
  missing_files_skullstripping = files_skullstripping[!file.exists(c(t1_fslbet_robust_path, brain_mask_path))]
  
  # only create directory if it doesn't already exist - EAH 5/5/26
  if (!dir.exists(brain.out.dir)) {
    dir.create(brain.out.dir,showWarnings = FALSE)
  }
  
  if ('t1_fslbet_robust' %in% names(missing_files_skullstripping)) {
      t1_fslbet_robust = fslbet_robust(t1_biascorrect,reorient = FALSE,correct = FALSE)
      writenii(t1_fslbet_robust, t1_fslbet_robust_path)
      writenii(t1_fslbet_robust, paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
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
reg.out.dir = paste0(main_path, "/data/", p, "/", s, "/registration/FLAIR_space")
if (argv$registration){
  # only create directory if it doesn't already exist - EAH 5/5/26
  if (!dir.exists(reg.out.dir)) {
    dir.create(reg.out.dir,showWarnings = FALSE, recursive = TRUE)
  }
  
  # file check - EAH 5/5/26
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
  t1_reg = readnii(paste0(reg.out.dir, "t1_n4_brain_reg_flair.nii.gz"))
  flair_n4_brain = readnii(paste0(reg.out.dir, "/flair_n4_brain.nii.gz"))
  brainmask_reg = readnii(paste0(reg.out.dir, "/brainmask_reg_flair"))
  if(!is.na(argv$t2)){
    t2_reg = readnii(paste0(main_path, "/data/", p, "/", s, "/registration/FLAIR_space/t2_n4_brain_reg_flair.nii.gz"))
  }
}

# WhiteStripe normalize data
white.out.dir = paste0(main_path, "/data/", p, "/", s, "/whitestripe/FLAIR_space")
if(argv$whitestripe){
  # only create directory if it doesn't already exist - EAH 5/5/26
  if (!dir.exists(white.out.dir)) {
    dir.create(white.out.dir,showWarnings = FALSE, recursive = TRUE)
  }
  
  # file check - EAH 5/5/26
  t1_n4_reg_brain_ws_path = paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws")
  t2_n4_reg_brain_ws_path = paste0(white.out.dir, "/t2_n4_brain_reg_flair_ws")
  flair_n4_brain_ws_path = paste0(white.out.dir, "/flair_n4_brain_ws")
  
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
} else{
  t1_n4_reg_brain_ws = readnii(paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws"))
  if(!is.na(argv$t2)){
    t2_n4_reg_brain_ws = readnii(paste0(white.out.dir, "/t2_n4_brain_reg_flair_ws"))
  }
  flair_n4_brain_ws = readnii(paste0(white.out.dir, "/flair_n4_brain_ws"))
}

# Mimosa
mim.out.dir = paste0(main_path, "/data/", p, "/", s, "/mimosa")

# only create directory if it doesn't already exist - EAH 5/5/26
if (!dir.exists(mim.out.dir)) {
  dir.create(mim.out.dir,showWarnings = FALSE)
}

# file check - EAH 5/5/26
probmap_path = paste0(mim.out.dir,"/mimosa_prob")
mimosa_mask_path = paste0(mim.out.dir,"/mimosa_mask")

files_mimosa = c(probmap = probmap_path, mimosa_mask = mimosa_mask_path)
missing_files_mimosa = files_mimosa[!file.exists(c(probmap_path, mimosa_mask_path))]

if ('probmap' %in% names(missing_files_mimosa)) {
  mimosa = mimosa_data(brain_mask=brainmask_reg, FLAIR=flair_n4_brain_ws, T1=t1_n4_reg_brain_ws, gold_standard=NULL, normalize="no", cores = 1, verbose = TRUE)
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
