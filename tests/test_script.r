# pacman::p_load(devtools)
# here::here()
# fname = here::here('../r/d/orig/S1_W1_D05_1.orig.bed.gz') 
# load_all()
# dt <- read_cpg_data(fname, "S1_W1_D05_1")

# pipeline <- process_pipeline("bedgraph")
# dt <- read_data(fname, pipeline)
# mdt <- 

# remove_leading_chr(dt)
# convert_meth_to_frac(dt, pipeline)
# collapse_strands = TRUE
# m_dt <- merge_strand_information(dt, pipeline, collapse_strands)
# if (! "cov" %in% names(m_dt)) {
#   m_dt[, cov := meth_cov + unmeth_cov]
#   m_dt[, .(meth_cov, unmeth_cov) := NULL]
# }
# upper_cov_cutoff = list(method = "percentile", value = 0.999)
# mask_technical_noise(m_dt, upper_cov_cutoff)
# reference = "hg19"
# align_to_reference = TRUE

# m_dt_all <- align_to_ref(m_dt, reference, align_to_reference)
# check()
