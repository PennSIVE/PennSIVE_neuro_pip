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


p <- arg_parser("Running lesion counting pipeline.", hide.opts = FALSE)
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
p <- add_argument(p, "--step", help = "Specify the step of lesion counting pipeline. preparation, count or consolidation.", default = "preparation")
p <- add_argument(p, "--method", help = "Specify the method of lesion counting. cc, jordan, or both.", default = "both")
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
    if(!is.na(argv$t2)){
      t2_biascorrect = bias_correct(file = t2,
                                 correction = "N4",
                                 verbose = TRUE)
      writenii(t2_biascorrect,paste0(bias.out.dir,"/T2_n4.nii.gz"))
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
    writenii(t1_fslbet_robust,paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
    }else{
      t1_fslbet_robust = readnii(paste0(bias.out.dir, "/T1_brain_n4.nii.gz"))
      } 
    }

  if (argv$skullstripping){
    dir.create(brain.out.dir,showWarnings = FALSE)
    t1_fslbet_robust = fslbet_robust(t1_biascorrect,reorient = FALSE,correct = FALSE)
    brain_mask = t1_fslbet_robust > 0 
    writenii(t1_fslbet_robust,paste0(brain.out.dir,"/T1_brain.nii.gz"))
    writenii(t1_fslbet_robust,paste0(bias.out.dir,"/T1_brain_n4.nii.gz"))
    writenii(brain_mask,paste0(brain.out.dir,"/T1_brainmask.nii.gz"))
  }

  # Registration to FLAIR Space
  reg.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/FLAIR_space")
  if (argv$registration){
    dir.create(reg.out.dir,showWarnings = FALSE, recursive = TRUE)
    ## Register T1 to FLAIR space 
    t1_to_flair = registration(filename = t1_biascorrect,
                             template.file = flair_biascorrect,
                             typeofTransform = "Rigid", remove.warp = FALSE,
                             outprefix=paste0(reg.out.dir,"/t1_reg_to_flair")) 

    t1_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(t1_fslbet_robust),
                                        transformlist = t1_to_flair$fwdtransforms, interpolator = "welchWindowedSinc"))
    brainmask_reg = ants2oro(antsApplyTransforms(fixed = oro2ants(flair_biascorrect), moving = oro2ants(brain_mask),
                                               transformlist = t1_to_flair$fwdtransforms, interpolator = "nearestNeighbor"))
    writenii(t1_reg, paste0(reg.out.dir,"/t1_n4_brain_reg_flair"))
    writenii(brainmask_reg, paste0(reg.out.dir,"/brainmask_reg_flair"))
    flair_n4_brain = flair_biascorrect
    flair_n4_brain[brainmask_reg==0] = 0
    writenii(flair_n4_brain, paste0(reg.out.dir,"/flair_n4_brain"))

    ## Register T2 to FLAIR space 
    if(!is.na(argv$t2)){
      t2_to_flair = registration(filename = t2_biascorrect,
                               template.file = flair_biascorrect,
                               typeofTransform = "Rigid", remove.warp = FALSE,
                               outprefix=paste0(reg.out.dir,"/t2_reg_to_flair"))
      t2_n4_brain = t2_to_flair$outfile
      t2_n4_brain[brainmask_reg==0] = 0
      writenii(t2_n4_brain, paste0(reg.out.dir,"/t2_n4_brain_reg_flair"))
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
    dir.create(white.out.dir,showWarnings = FALSE, recursive = TRUE)
    ind1 = whitestripe(t1_reg, "T1")
    t1_n4_reg_brain_ws = whitestripe_norm(t1_reg, ind1$whitestripe.ind)
    writenii(t1_n4_reg_brain_ws, paste0(white.out.dir,"/t1_n4_brain_reg_flair_ws"))
    if(!is.na(argv$t2)){
      ind2 = whitestripe(t2_n4_brain, "T2")
      t2_n4_reg_brain_ws = whitestripe_norm(t2_n4_brain, ind2$whitestripe.ind)
      writenii(t2_n4_reg_brain_ws, paste0(white.out.dir,"/t2_n4_brain_reg_flair_ws"))
    }
    ind3 = whitestripe(flair_n4_brain, "T2")
    flair_n4_brain_ws = whitestripe_norm(flair_n4_brain, ind3$whitestripe.ind)
    writenii(flair_n4_brain_ws, paste0(white.out.dir,"/flair_n4_brain_ws"))
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
    writenii(probmap > as.numeric(argv$threshold), paste0(mim.out.dir,"/mimosa_mask"))
  }else{
    probmap = readnii(paste0(mim.out.dir,"/mimosa_prob"))
  }

  # Register to EPI Space
  reg.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/EPI_space")
  if(!file.exists(reg.epi.out.dir)){
    dir.create(reg.epi.out.dir, showWarnings = FALSE)
    # Fast
    t1_fast = fast(t1_reg, bias_correct = F, opts = "-t 1 -n 3")
    flair_to_epi = registration(filename = flair_biascorrect,
                                  template.file = abs(phase_biascorrect),
                                  typeofTransform = "Rigid", remove.warp = FALSE) ### rigid

    brainmask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_biascorrect)), moving = oro2ants(brainmask_reg),
                                                    transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
    writenii(brainmask_reg_epi, paste0(reg.epi.out.dir,'/brainmask_reg_epi'))
    phase_n4_brain = phase_biascorrect * brainmask_reg_epi
    writenii(phase_n4_brain, paste0(reg.epi.out.dir,'/phase_n4_brain'))

    t1_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(t1_reg),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
    flair_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(flair_n4_brain),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
    mimosa_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(probmap),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "welchWindowedSinc"))
    mimosa_mask_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(probmap>0.2),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))
    fast_reg_epi = ants2oro(antsApplyTransforms(fixed = oro2ants(abs(phase_n4_brain)), moving = oro2ants(t1_fast),
                transformlist = flair_to_epi$fwdtransforms, interpolator = "nearestNeighbor"))

    writenii(t1_reg_epi, paste0(reg.epi.out.dir, "/t1_reg_epi"))
    writenii(flair_reg_epi, paste0(reg.epi.out.dir, "/flair_reg_epi"))
    writenii(mimosa_reg_epi, paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
    writenii(mimosa_mask_reg_epi, paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
    writenii(fast_reg_epi, paste0(reg.epi.out.dir,'/fast_reg_epi'))
  }else{
    mimosa_reg_epi = readnii(paste0(reg.epi.out.dir, "/mimosa_reg_epi"))
    mimosa_mask_reg_epi = readnii(paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
    fast_reg_epi = readnii(paste0(reg.epi.out.dir,'/fast_reg_epi'))
  }

  if(argv$method == 'jordan' | argv$method == 'both'){
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
  }


}else if(argv$step == "count"){
  p = argv$participant
  ses = argv$session
  reg.epi.out.dir = paste0(main_path, "/data/", p, "/", ses, "/registration/EPI_space")
  count.out.dir = paste0(main_path, "/data/", p, "/", ses, "/count")
  dir.create(count.out.dir,showWarnings = FALSE)
    if(argv$method == 'cc'){
        mimosa_mask_reg_epi = readnii(paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
        lesion_num_cc = max(ants2oro(labelClusters(oro2ants(mimosa_mask_reg_epi > 0),minClusterSize=27)))
        write.csv(lesion_num_cc,paste0(count.out.dir,"/",p,"_connected_components.csv"))        
    
    }else if(argv$method == 'jordan'){
        label_result = readnii(paste0(reg.epi.out.dir,"/lesions_reg_epi_labeled"))
        lesion_num_jordan = max(label_result)
        write.csv(lesion_num_jordan,paste0(count.out.dir,"/",p,"_jordan.csv"))
    
    }else if(argv$method == 'both'){
        mimosa_mask_reg_epi = readnii(paste0(reg.epi.out.dir, "/mimosa_mask_reg_epi"))
        lesion_num_cc = max(ants2oro(labelClusters(oro2ants(mimosa_mask_reg_epi > 0),minClusterSize=27)))
        write.csv(lesion_num_cc,paste0(count.out.dir,"/",p,"_connected_components.csv"))
        
        label_result = readnii(paste0(reg.epi.out.dir,"/lesions_reg_epi_labeled"))
        lesion_num_jordan = max(label_result)
        write.csv(lesion_num_jordan,paste0(count.out.dir,"/",p,"_jordan.csv"))
    }


}else if(argv$step == "consolidation"){
    if(argv$method == 'cc' | argv$method == 'both'){
        cc_files = list.files(paste0(main_path, "/data"), pattern = "*_connected_components.csv", recursive = TRUE, full.names = TRUE) 
        cc_con = lapply(cc_files, function(x){
            sub_file = read_csv(x)
            subj = str_split(x, "/")[[1]][length(str_split(x, "/")[[1]]) - 3]
            sub_file = sub_file %>% mutate(subject = subj)
        }) %>% bind_rows()
        cc_con = cc_con[,-1]
        colnames(cc_con) = c("lesion_count", "subject")
        if(!file.exists(paste0(main_path, "/stats"))){
            dir.create(paste0(main_path, "/stats"))
        }
        write_csv(cc_con, paste0(main_path, "/stats/lesion_count_connected_components.csv"))
    }
    if(argv$method == 'jordan' | argv$method == 'both') {
        jordan_files = list.files(paste0(main_path, "/data"), pattern = "*_jordan.csv", recursive = TRUE, full.names = TRUE) 
        jordan_con = lapply(jordan_files, function(x){
            sub_file = read_csv(x)
            subj = str_split(x, "/")[[1]][length(str_split(x, "/")[[1]]) - 3]
            sub_file = sub_file %>% mutate(subject = subj)
        }) %>% bind_rows()
        jordan_con = jordan_con[,-1]
        colnames(jordan_con) = c("lesion_count", "subject")
        if(!file.exists(paste0(main_path, "/stats"))){
            dir.create(paste0(main_path, "/stats"))
        }
        write_csv(jordan_con, paste0(main_path, "/stats/lesion_count_jordan.csv"))
    }
}



# }else if(argv$step == "consolidation"){
#   prl_files = list.files(paste0(main_path, "/data"), pattern = "*_preds.csv", recursive = TRUE, full.names = TRUE) 
#   prl_con = lapply(prl_files, function(x){
#     sub_file = read_csv(x)
#     subj = str_split(x, "/")[[1]][length(str_split(x, "/")[[1]]) - 2]
#     sub_file = sub_file %>% mutate(subject = subj)
#   }) %>% bind_rows()
#   colnames(prl_con) = c("lesion_id", "rim_neg", "rim_pos", "subject")
#   if(!file.exists(paste0(main_path, "/stats"))){
#     dir.create(paste0(main_path, "/stats"))
#   }
#   write_csv(prl_con, paste0(main_path, "/stats/prl_probability.csv"))
# }
