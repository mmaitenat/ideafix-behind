# This code extracts info from a vcf file. 
ExtractFromVcf <- function(vcf.filename, 
                           sample.id,
                           read.length,
                           is.paired,
                           allele.freq = NULL, 
                           alt.bases = NULL, 
                           ref.bases = NULL, 
                           ref.allele = NULL, 
                           alt.allele = NULL,
                           base.qual = NULL,
                           base.qual.frac = NULL,
                           frag.length = NULL,
                           pos.from.end = NULL,
                           map.qual = NULL, 
                           FdeamC = NULL,
                           SB = NULL) {
  
  # Check if any of the arguments is TRUE. Exit if not.
  arguments <- as.list(environment())
  if(all(sapply(arguments[-(1:4)], is.null))) {
    stop ("At least one vcf feature should be asked for.")
  }
  
  # Define functions for each argument  
  args.to.vcf.funs <- list(
    allele.freq = "GetAF",
    alt.bases = "GetALTBases",
    ref.bases = "GetREFBases",
    ref.allele = "GetREFAllele",
    alt.allele = "GetALTAllele",
    base.qual = "GetBaseQual",
    base.qual.frac = "GetBaseQualFrac",
    frag.length = "GetFragLength",
    pos.from.end = "GetPosFromEnd",
    map.qual = "GetMapQual",
    FdeamC = "GetFdeamC",
    SB = "GetSB")
  
  # Obtain data from vcf
  prov.arguments <- arguments[!sapply(arguments, is.null)]
  ## Field. Depending on if it's a paired variant calling or a single one, the field we are interested in is the 11th or 10th
  # field <- ifelse(is.paired, 11, 10)
  ## Depth for each position. We calculate this in the outside because it's used more than once and this way we optimize the calculations
  depth <- GetDepth(vcf.filename, sample.id)
  out.data <- sapply(names(prov.arguments[-(1:4)]), function(x) {
    fun <- get(args.to.vcf.funs[[x]])
    fun(vcf.filename, field = field, depth = depth, read.length = read.length, sample.id = sample.id)
  }
  )
  out.data <- c(out.data, id = list(GetMutID(vcf.filename)))
  out.data <- lapply(rapply(out.data, enquote, how = "unlist"), eval)
  X <- bind_rows(out.data)
  # See if variant is a SNP and, based on that, create a binary column (isSNP) that will only be 1 if variant is present in dbSNP
  X %>%
    mutate(isSNP = isSNP(vcf.filename)) -> X
  # Add read-length as a feature:
  X %>%
    mutate(read.length = read.length) -> X
  # Create binary feature "both reads aligned":
  X %>%
    mutate(both.reads.aligned = ifelse(frag.length1 != 0, 1, 0)) -> X
  # Add the variable "normalized distance from end" by combining two columns
  if (!is.null(pos.from.end)) { # if we have asked for that variable
    X %>%
      mutate(norm.pos.from.end = pos.from.end/frag.length1) -> X
  }
  # Rename columns
  rename(X, norm.alt.bases = alt.bases.total.depth,
         norm.ref.bases = ref.bases.total.depth,
         frag.length = frag.length1,
         frag.length.frac = frag.length2,
         FdeamC = FdeamC1,
         SOB = FdeamC2,
         SBGuo = SB1,
         SBGATK = SB2) -> X
}

GetDepth <- function(vcf.filename, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD]\n'", sample.id, vcf.filename) # FORMAT fields need the brackets, INFO dont
  depth <- system(cmd, intern = TRUE) %>%
    as_tibble() %>%
    separate(value, into = c("ref", "alt"), sep = ",") %>% # this will produce NAs when we have tri/cuatriallelic sites. However it's not a problem because we later discard those sites and only keep C:G>T:A
    mutate_all(funs(as.numeric)) %>%
    mutate(total.depth = ref + alt) %>%
    select(total.depth)
}

GetMutID <- function(vcf.filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%CHROM:%%POS\n' %s", vcf.filename)
  id <- system(cmd, intern = TRUE)
  return(id)
}

GetAF <- function(vcf.filename, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AF]\n'", sample.id, vcf.filename)
  AF <- system(cmd, intern = TRUE)
  AF <- as.numeric(AF)
  return(AF)
}

GetREFAllele <- function(vcf.filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%REF\n' %s", vcf.filename)
  REF.alleles <- system(cmd, intern = TRUE)
  return(REF.alleles)
}

GetALTAllele <- function(vcf.filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%ALT\n' %s", vcf.filename)
  ALT.alleles <- system(cmd, intern = TRUE)
  return(ALT.alleles)
}

GetBaseQual <- function(vcf.filename, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MBQ{1}]\n'", sample.id, vcf.filename)
  base.qual <- system(cmd, intern = TRUE)
  base.qual <- as.numeric(base.qual)
  return(base.qual)
}

GetBaseQualFrac <- function(vcf.filename, sample.id,  ...) {
  cmd.tumor <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MBQ{1}]\n'", sample.id, vcf.filename)
  cmd.normal <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MBQ{0}]\n'", sample.id, vcf.filename)
  base.qual.tumor <- system(cmd.tumor, intern = TRUE) %>%
    as.numeric()
  base.qual.normal <- system(cmd.normal, intern = TRUE) %>%
    as.numeric()
  base.qual.frac <- base.qual.tumor/base.qual.normal
  return(base.qual.frac)
}

GetFragLength <- function(vcf.filename, read.length, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MFRL{1}]\n'", sample.id, vcf.filename)
  # Median fragment length for reads that support the alternate allele
  fraglen <- system(cmd, intern = TRUE) %>%
    as.numeric()
  fraglen.ratio <- fraglen/read.length
  return(list(fraglen, fraglen.ratio))
}

GetPosFromEnd <- function(vcf.filename, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MPOS]\n'", sample.id, vcf.filename)
  dist.from.end <- system(cmd, intern = TRUE) %>%
    as.numeric() # NAs will appear for multiallelic variants
  return(dist.from.end)
}

GetFdeamC <- function(vcf.filename, sample.id, ...) {
  ref.allele <- GetREFAllele(vcf.filename)
  alt.allele <- GetALTAllele(vcf.filename)
  F1R2.alt.cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%F1R2{1}]\n'", sample.id, vcf.filename)
  F1R2.alt <- system(F1R2.alt.cmd, intern = TRUE) %>%
    as.numeric()
  F2R1.alt.cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%F2R1{1}]\n'", sample.id, vcf.filename)
  F2R1.alt <- system(F2R1.alt.cmd, intern = TRUE) %>%
    as.numeric()
  denom <- F2R1.alt + F1R2.alt
  numer <- F1R2.alt # default numerator
  GA.idx <- ref.allele == "G" & alt.allele == "A"
  numer[GA.idx] <- F2R1.alt[GA.idx] # change if mutation is G>A
  FdeamC <- numer/denom
  # does not apply to not FFPE mutations
  notFFPE.idx <- !(paste(ref.allele, alt.allele, sep = ":") %in% c("C:T", "G:A"))
  FdeamC[notFFPE.idx] <- NA
  SOB <- (F1R2.alt - F2R1.alt) /denom
  return(list(FdeamC, SOB))
}

GetMapQual <- function(vcf.filename, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MMQ]\n'", sample.id, vcf.filename)
  map.qual <- system(cmd, intern = TRUE) %>%
    as.numeric() # NAs will appear for multiallelic variants
  return(map.qual)
}

GetALTBases <- function(vcf.filename, depth, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD{1}]\n'", sample.id, vcf.filename)
  ALT.bases <- system(cmd, intern = TRUE)
  ALT.bases <- as.numeric(ALT.bases)
  ALT.base.ratio <- ALT.bases/depth
  return(list(ALT.bases, ALT.base.ratio))
}

GetREFBases <- function(vcf.filename, depth, sample.id, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD{0}]\n'", sample.id, vcf.filename)
  REF.bases <- system(cmd, intern = TRUE)
  REF.bases <- as.numeric(REF.bases)
  REF.base.ratio <- REF.bases/depth
  return(list(REF.bases, REF.base.ratio))
}

GetSB <- function(vcf.filename, sample.id, ...) {
  # This function will provide 2 values: 
  # FS: phred-scaled p-value using Fisher's exact test to detect strand bias (NOT ANYMORE)
  # SB: SB as defined in  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/ . Ranges from 0 to infinity.
  # GATKSB: SB as defined by GATK and described in  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/ (same as above). Ranges from 0 to infinity.
  cmd.SB <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%SB]\n'", sample.id, vcf.filename)
 # cmd.SOR <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%INFO/SOR\n' %s", vcf.filename)
 # cmd.FS <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%INFO/FS\n' %s", vcf.filename)
 # FS <- system(cmd.FS, intern = TRUE) %>% 
 # as.numeric()
 # SOR <- system(cmd.SOR, intern = TRUE) %>%
 # as.numeric()
 SB.all <- system(cmd.SB, intern = TRUE) %>%
   as_tibble() %>%
   separate(value, into = c("REF_fw", "REF_rev", "ALT_fw", "ALT_rev"), sep = ",") %>%
   mutate_all(funs(as.numeric)) %>%
   mutate_all(funs(replace(., . == 0, 0.0001))) %>% # I am putting pseudozeros not to divide by 0
   mutate(SB = CalcSBGuo(a = REF_fw, b = ALT_fw, e = REF_rev, d = ALT_rev),
          SBGATK = CalcSBGATK(a = REF_fw, b = ALT_fw, e = REF_rev, d = ALT_rev))
 return(list(pull(SB.all, SB), pull(SB.all, SBGATK)))
}

isSNP <- function(vcf.filename, ...) {
  cmd <- sprintf("cat %s | awk '$7==\"PASS\"' | awk '{if ($3 != \".\") print 1; else print 0;}'", vcf.filename)
  SNP.id <- system(cmd, intern = TRUE) %>%
    as.factor()
  return(SNP.id)
}

# The following two definitions for strand bias were obtained from here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/
CalcSBGuo <- function(a, b, e, d) { # Guo because it's the first author of the article to which the previous one references this calculation
  denom1 <- a + b
  denom2 <- e + d
  denom3 <- a + b + e + d
  numer <- abs((b / denom1) - (d / denom2))
  denom <- (b + d) / denom3
  return(numer / denom)
}

CalcSBGATK <- function (a, b, e, d) {
  numer1 <- (b * e) / ((a + b) * (e + d))
  numer2 <- (d * a) / ((e + d) * (a + b))
  denom <- (a + e) / (a + b + e + d)
  return(pmax(numer1/denom, numer2/denom))
}
