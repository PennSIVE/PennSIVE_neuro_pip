suppressMessages(library(argparser))
suppressMessages(library(neurobase))
suppressMessages(library(ALPaCA))
suppressMessages(library(extrantsr))
suppressMessages(library(WhiteStripe))
suppressMessages(library(fslr))
suppressMessages(library(ANTsR))
suppressMessages(library(ANTsRCore))
suppressMessages(library(mimosa))
suppressMessages(library(torch))
suppressMessages(library(purrr))

p <- arg_parser("Running ALPaCA pipeline to assess the CVS, PRLs, and lesions.", hide.opts = FALSE)
p <- add_argument(p, "--mainpath", short = '-m', help = "Specify the main path where MRI images can be found.")
p <- add_argument(p, "--participant", short = '-p', help = "Specify the subject id.")
p <- add_argument(p, "--session", short = '-ses', help = "Specify the session id.")
p <- add_argument(p, "--t1", help = "Specify the T1 sequence name.")
p <- add_argument(p, "--flair", help = "Specify the FLAIR sequence name.")
p <- add_argument(p, "--epi_mag", help = "Specify the EPI sequence magnitude image name.")
p <- add_argument(p, "--epi_phase", help = "Specify the EPI sequence unwrapped phase image name.")
p <- add_argument(p, "--n4", help = "Specify whether to run bias correction step.", default = TRUE)
p <- add_argument(p, "--skullstripping", short = '-s', help = "Specify whether to run skull stripping step.", default = TRUE)
p <- add_argument(p, "--registration", short = '-r', help = "Specify whether to run registration to FLAIR space step.", default = TRUE)
p <- add_argument(p, "--whitestripe", short = '-w', help = "Specify whether to run whitestripe step.", default = TRUE)
p <- add_argument(p, "--mimosa", help = "Specify whether to run mimosa segmentation step.", default = TRUE)
p <- add_argument(p, "--threshold", help = "Specify the threshold used to generate mimosa candidate mask.", default = 0.05)
p <- add_argument(p, "--step", help = "Specify the step of cvs pipeline. estimation or consolidation.", default = "estimation")
p <- add_argument(p, "--lesioncenter", help = "Provide the path to the lesioncenter package.")
p <- add_argument(p, "--hdbetpath", help = "Specify the path to the HD-BET binary", default = "~/.local/bin/hd-bet")
p <- add_argument(p, "--mpath", help = "Specify the path to the trained mimosa model.")
p <- add_argument(p, "--helpfunc", help = "Specify the path to the help functions.")
argv <- parse_args(p)

# Read in Files
main_path = argv$mainpath

if (argv$step == "estimation"){
  # Load lesion center package
  my_path = paste0(argv$lesioncenter, "/")
  source_files = list.files(my_path)
  source(paste0(argv$helpfunc, "/label_code.R"))

  # Read and check inputs
  purrr::map(paste0(my_path, source_files), source)
  p = argv$participant
  ses = argv$session
  anat.path<-paste0(main_path, "/data/", p,  "/", ses, "/anat/")
  message('Checking inputs...')
  if(is.na(argv$t1)) stop("Missing T1 sequence!")else{
    t1 = readnii(paste0(anat.path, argv$t1))
    }
    if(is.na(argv$flair)) stop("Missing FLAIR sequence!")else{
      flair = readnii(paste0(anat.path, argv$flair))
    }

    if(is.na(argv$epi_mag)) stop("Missing EPI magnitude image!")else{
      epi_map = read_rpi(paste0(anat.path, argv$epi_mag))
    }
    if(is.na(argv$epi_phase)) stop("Missing EPI phase image!")else{
      epi_phase = read_rpi(paste0(anat.path, argv$epi_phase))
    }

  # Specify output directory
  outdir = paste0(main_path, "/data/", p,  "/", ses, "/alpaca")
  dir.create(outdir,showWarnings = FALSE)

  ## Function to run HD-BET from R

  hdbet_path = argv$hdbetpath
  hdbet <- function(img,
                    mask = NULL,
                    device = NULL,
                    cleanup = TRUE,
                    verbose = FALSE,
                    hdbet_bin = hdbet_path) {
    hdbet_bin <- path.expand(hdbet_bin)
    img_file <- neurobase::checkimg(img)
    if (!file.exists(hdbet_bin)) {
      stop(
        "Cannot find hd-bet executable at: ", hdbet_bin, "\n",
        "Either install hd-bet there, or pass `hdbet_bin` with the correct path."
      )
    }
    out_base <- tempfile()
    out_file <- paste0(out_base, ".nii.gz")
    mask_base <- if (!is.null(mask)) tempfile() else NULL
    mask_file <- if (!is.null(mask)) paste0(mask_base, ".nii.gz") else NULL
    cmd <- paste(
      shQuote(hdbet_bin),
      "-i", shQuote(img_file),
      "-o", shQuote(out_file),
      if (!is.null(mask)) paste("-m", shQuote(mask_file)) else "",
      if (!is.null(device)) paste("-device", shQuote(device),"--disable_tta") else ""
    )
    if (verbose) message(cmd)
    status <- system(cmd)
    if (status != 0) stop("HD-BET failed (exit status ", status, ")")
    brain <- oro.nifti::readNIfTI(out_base, reorient = FALSE)
    if (is.null(mask)) {
      if (cleanup) unlink(out_file)
      return(brain)
    }
    brain_mask <- oro.nifti::readNIfTI(mask_base, reorient = FALSE)
    if (cleanup) unlink(c(out_file, mask_file))
    list(brain = brain, mask = brain_mask)
  }


  message('Loading inputs...')
  t1_path = paste0(anat.path, argv$t1)
  flair_path = paste0(anat.path, argv$flair)
  epi_path = paste0(anat.path, argv$epi_mag)
  phase_path = paste0(anat.path, argv$epi_phase)

  # Read in images
  t1 <- oro2ants(read_rpi(t1_path, verbose = verbose))
  flair <- oro2ants(read_rpi(flair_path, verbose = verbose))
  epi <- oro2ants(read_rpi(epi_path, verbose = verbose))
  phase <- oro2ants(read_rpi(phase_path, verbose = verbose))

  # Bias Correction
  message('Bias correction...')
  bias.out.dir = paste0(main_path, "/data/", p, "/", ses, "/bias_correction")
  if(argv$n4){
    dir.create(bias.out.dir,showWarnings = FALSE)
    flair_biascorrect = bias_correct(file = flair,
                                      correction = "N4",
                                      verbose = TRUE)
    writenii(flair_biascorrect,paste0(bias.out.dir,"/FLAIR_n4.nii.gz"))
    t1_biascorrect = bias_correct(file = t1,
                                  correction = "N4",
                                  verbose = TRUE)
    writenii(t1_biascorrect,paste0(bias.out.dir,"/T1_n4.nii.gz"))
    phase_biascorrect = bias_correct(file = phase,
                                  correction = "N4",
                                  verbose = TRUE)
    writenii(phase_biascorrect,paste0(bias.out.dir,"/PHASE_n4.nii.gz"))
    epi_biascorrect = bias_correct(file = epi,
                                  correction = "N4",
                                  verbose = TRUE)
    writenii(epi_biascorrect,paste0(bias.out.dir,"/EPI_n4.nii.gz"))
      
  }else{
    t1_biascorrect = readnii(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/T1_n4.nii.gz"))
    flair_biascorrect = readnii(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/FLAIR_n4.nii.gz"))
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
    epi_file = existing_files[which(grepl("EPI_n4*", existing_files))]
    if(length(epi_file) == 0){
      epi_biascorrect = bias_correct(file = epi,
                                  correction = "N4",
                                  verbose = TRUE)
      writenii(epi_biascorrect,paste0(bias.out.dir,"/EPI_n4.nii.gz"))
    }else{
      epi_biascorrect = read_rpi(paste0(main_path, "/data/", p, "/", ses, "/bias_correction/EPI_n4.nii.gz"))
    }
  }

  # Skull Stripping
  message('Skull stripping...')
  brain.out.dir = paste0(main_path, "/data/", p,  "/", ses, "/t1_brain")
  if(!argv$skullstripping){
    brain_paths = list.files(brain.out.dir, recursive = TRUE, full.names = TRUE)
    brain_mask_path = brain_paths[which(grepl("*brainmask.nii.gz$", brain_paths))]
    brain_mask = readnii(brain_mask_path)
    t1_brain = t1_biascorrect * brain_mask
    writenii(t1_brain, paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
  }

  if (argv$skullstripping){
    dir.create(brain.out.dir,showWarnings = FALSE, recursive = TRUE)
    message("Running HD-BET...")
    t1_brain<-hdbet(t1_biascorrect,device="cpu",hdbet_bin = hdbet_path)
    message("Finished HD-BET...")
    brain_mask = t1_brain > 0 
    writenii(t1_brain,paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
    writenii(brain_mask,paste0(brain.out.dir,"/T1_brainmask.nii.gz"))
  }


  # Registration to FLAIR space
  message('FLAIR space registration...')
  reg.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/FLAIR_space")
  if (argv$registration){
      dir.create(reg.out.dir,showWarnings = FALSE, recursive = TRUE)
      ## Register T1 to FLAIR space 
      t1_to_flair = registration(filename = t1_biascorrect,
                              template.file = flair_biascorrect,
                              typeofTransform = "Rigid", remove.warp = FALSE,
                              outprefix=paste0(reg.out.dir,"/t1_reg_to_flair")) 

      t1_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(t1_brain),
                                          transformlist = t1_to_flair$fwdtransforms, interpolator = "welchWindowedSinc"))
      brainmask_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(brain_mask),
                                                transformlist = t1_to_flair$fwdtransforms, interpolator = "nearestNeighbor"))
      writenii(t1_reg, paste0(reg.out.dir,"/t1_n4_brain_reg_flair"))
      writenii(brainmask_reg, paste0(reg.out.dir,"/brainmask_reg_flair"))
      flair_n4_brain = flair_biascorrect
      flair_n4_brain[brainmask_reg==0] = 0
      writenii(flair_n4_brain, paste0(reg.out.dir,"/flair_n4_brain"))
  }else{
      t1_reg = readnii(paste0(reg.out.dir, "/t1_n4_brain_reg_flair.nii.gz"))
      flair_n4_brain = readnii(paste0(reg.out.dir, "/flair_n4_brain.nii.gz"))
      brainmask_reg = readnii(paste0(reg.out.dir, "/brainmask_reg_flair"))
  }

  # WhiteStripe normalize data
  message('WhiteStripe normalization...')
  white.out.dir = paste0(main_path, "/data/", p, "/", ses, "/whitestripe/FLAIR_space")
  if(argv$whitestripe){
      dir.create(white.out.dir,showWarnings = FALSE, recursive = TRUE)
      ind1 = whitestripe(t1_reg, "T1")
      t1_n4_reg_brain_ws = whitestripe_norm(t1_reg, ind1$whitestripe.ind)
      writenii(t1_n4_reg_brain_ws, paste0(white.out.dir,"/t1_n4_brain_reg_flair_ws"))
      ind3 = whitestripe(flair_n4_brain, "T2")
      flair_n4_brain_ws = whitestripe_norm(flair_n4_brain, ind3$whitestripe.ind)
      writenii(flair_n4_brain_ws, paste0(white.out.dir,"/flair_n4_brain_ws"))
  }else{
        t1_n4_reg_brain_ws = readnii(paste0(white.out.dir, "/t1_n4_brain_reg_flair_ws"))
        flair_n4_brain_ws = readnii(paste0(white.out.dir, "/flair_n4_brain_ws"))
  }

  # Mimosa
  message('MIMoSA segmentation...')
  mim.out.dir = paste0(main_path, "/data/", p, "/", ses, "/mimosa")
  if(argv$mimosa){
      dir.create(mim.out.dir,showWarnings = FALSE)

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
      writenii(probmap > as.numeric(argv$threshold), paste0(mim.out.dir,"/mimosa_cand_mask"))
  }else{
      probmap = readnii(paste0(mim.out.dir,"/mimosa_prob"))
  }

  # Threshold MIMoSA mask and identify/split confluent lesions
  prob_05 <- antsImageClone(oro2ants(probmap) > as.numeric(argv$threshold))
  if (sum(prob_05) == 0) {
      prob_05_labeled <- antsImageClone(prob_05)
      prob_05_erode <- antsImageClone(prob_05)
  } else {
      prob_05_labeled <- oro2ants(ALPaCA:::label_lesion(probmap, prob_05, mincluster = 30))
      prob_05_erode <- iMath(prob_05_labeled, "GE", 1)
  }
  antsImageWrite(prob_05_labeled, file.path(mim.out.dir, "/labeled_candidates_flairspace.nii.gz"))
  antsImageWrite(prob_05_erode, file.path(mim.out.dir, "/eroded_candidates_flairspace.nii.gz"))

  # Register to EPI Space
  message('Registration to EPI space...')
  reg.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/EPI_space")

  reg.epi.files <- paste0(reg.epi.out.dir, c(
    "/t1_reg_epi.nii.gz",
    "/flair_reg_epi.nii.gz",
    "/mimosa_reg_epi.nii.gz",
    "/mimosa_cand_mask_reg_epi.nii.gz",
    "/labeled_candidates.nii.gz",
    "/eroded_candidates.nii.gz",
    "/epi_n4_brain",
    "phase_n4_brain",
    "brainmask_reg_epi"
  ))

  if(!all(file.exists(reg.epi.files))){
      dir.create(reg.epi.out.dir, showWarnings = FALSE)
      flair_to_epi = registration(filename = flair_biascorrect,
                                    template.file = abs(phase_biascorrect),
                                    typeofTransform = "Rigid", remove.warp = FALSE) ### rigid

      brainmask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_biascorrect)), moving = oro2ants(brainmask_reg),
                                                      transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      writenii(brainmask_reg_epi, paste0(reg.epi.out.dir,'/brainmask_reg_epi'))
      phase_n4_brain = phase_biascorrect * brainmask_reg_epi
      writenii(phase_n4_brain, paste0(reg.epi.out.dir,'/phase_n4_brain'))
      epi_n4_brain = epi_biascorrect * brainmask_reg_epi
      writenii(epi_n4_brain, paste0(reg.epi.out.dir,'/epi_n4_brain'))

      t1_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(t1_reg),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
      flair_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(flair_n4_brain),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
      mimosa_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(probmap),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
      mimosa_mask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(probmap>as.numeric(argv$threshold)),
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      labeled_candidates = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = prob_05_labeled,
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      eroded_candidates = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = prob_05_erode,
                  transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
      
      writenii(t1_reg_epi, paste0(reg.epi.out.dir, "/t1_reg_epi"))
      writenii(flair_reg_epi, paste0(reg.epi.out.dir, "/flair_reg_epi"))
      writenii(mimosa_reg_epi, paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
      writenii(mimosa_mask_reg_epi, paste0(reg.epi.out.dir, "/mimosa_cand_mask_reg_epi"))
      writenii(labeled_candidates, paste0(reg.epi.out.dir, "/labeled_candidates"))
      writenii(eroded_candidates, paste0(reg.epi.out.dir, "/eroded_candidates"))
  }else{
      t1_reg_epi <- readnii(paste0(reg.epi.out.dir, "/t1_reg_epi"))
      flair_reg_epi <- readnii(paste0(reg.epi.out.dir, "/flair_reg_epi"))
      mimosa_reg_epi <- readnii(paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
      mimosa_mask_reg_epi <- readnii(paste0(reg.epi.out.dir, "/mimosa_cand_mask_reg_epi"))
      labeled_candidates <- readnii(paste0(reg.epi.out.dir, "/labeled_candidates"))
      eroded_candidates <- readnii(paste0(reg.epi.out.dir, "/eroded_candidates"))
  }

  message('Making predictions...')
  make_predictions(
  t1 = t1_reg_epi,
  flair = flair_reg_epi,
  epi = epi_n4_brain,
  phase = phase_n4_brain,
  labeled_candidates = labeled_candidates,
  eroded_candidates = eroded_candidates,
  output_dir = outdir
  )
} else if(argv$step == "consolidation"){
  alpaca_con = list.files(paste0(main_path, "/data"), pattern = "probabilities.csv", recursive = TRUE, full.names = TRUE) %>% read_csv() %>% bind_rows()
  if(!file.exists(paste0(main_path, "/stats"))){
    dir.create(paste0(main_path, "/stats"))
  }
  write_csv(alpaca_con, paste0(main_path, "/stats/alpaca_score.csv"))
}