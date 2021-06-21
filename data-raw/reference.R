## code to prepare `reference/hg37` dataset goes here


dt1 <- get_ref_cpgs("hg19", "BSgenome.Hsapiens.UCSC.hg19",
  dir = fs::path("inst", "extdata"))
dt2 <- get_ref_cpgs("hg38", "BSgenome.Hsapiens.UCSC.hg38",
  dir = fs::path("inst", "extdata"))


dt1 <- get_ref_cpgs("hg19", "BSgenome.Hsapiens.UCSC.hg19")
dt2 <- get_ref_cpgs("hg38", "BSgenome.Hsapiens.UCSC.hg38")
