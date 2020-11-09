ExtractFromFa <- function(vcf.filename = NULL,
                          mut.id = NULL,
                          fa.filename, 
                          k = 1) {
  
  # vcf filename of mutation id list should be provided
  if (is.null(vcf.filename) & is.null(mut.id)) {
    stop("vcf filename of mutation id list should be provided.")
  }

  sprintf("Variant context of +- %s positions will be extracted.", k)
  
  # If id list is not provided, obtain it from vcf
  if (is.null(mut.id)) {
    mut.id <- GetMutID(vcf.filename)
  }
  
  surr.bases <- GetSurrBases(fa.filename, mut.id, k)
  # Formatting
  odd.idx <- seq(from = 1, to = length(surr.bases), by = 2)
  even.idx <- seq(from = 2, to = length(surr.bases), by = 2)
  base.datamat <- data.frame(region = surr.bases[odd.idx], 
                             bases = surr.bases[even.idx]) %>%
    as_tibble()
  # Separate bases in before and after & obtain mutation id_from region
  base.datamat %>%
    separate(col = bases, into = c("before", "after"), sep = 3, remove = FALSE) %>%
    separate(col = before, into = c("before", "current"), sep = -1, remove = TRUE) %>%
    separate(col = before, into = c("before.2", "before.1"), sep = 1, remove = FALSE) %>%
    separate(col = after, into = c("after.1", "after.2"), sep = 1, remove = FALSE) %>%
    mutate(chr = sub(":.*", "", region) %>%
             sub(">", "", .)) %>%
    mutate(pos = sub(".*:", "", region) %>% 
             sub("-.*", "", .) %>%
             as.numeric() %>%
             magrittr::add(k)) %>%
    mutate(id = paste(chr, pos, sep = ":")) %>%
    select(id, current, bases, before.2, before.1, after.1, after.2, before, after) -> base.datamat
   # Obtain variables is.repeat.region. The variable will be TRUE if any of the position before, the position itself or the position after are lowercase (= repeat region)
  base.datamat %>%
    mutate(is.repeat.region = base.datamat %>%
             select(c("current", "before.1", "after.1")) %>%
             pmap_lgl(~any(checkLowerCase(.)))) -> base.datamat
  
  # Convert to uppercase
  cols.toupper <- setdiff(colnames(base.datamat), c("id", "bases", "is.repeat.region"))
  base.datamat %>%
    mutate_at(cols.toupper, toupper) -> base.datamat
  base.datamat
}

GetSurrBases <- function(fa.filename, mut.id, k) {
  tmp.dir <- "/tmp"
  tmp.filename <- file.path(tmp.dir, "mut_regions.temp")
  WriteRegionFile(mut.id, k, tmp.filename)
  cmd <- sprintf("samtools faidx %s -r %s", fa.filename, tmp.filename)
  bases <- system(cmd, intern = TRUE)
  return(bases)
}

# Helper functions
WriteRegionFile <- function(mut.list, k, tmp.filename) {
  strsplit(mut.list, ":") %>% 
    data.frame() %>%
    t() %>%
    data.frame() %>%
    magrittr::set_colnames(c("chr", "pos")) %>%
    mutate_at(1, as.character) %>%
    mutate_at(2, as.character) %>%
    mutate_at(2, as.numeric) %>%
    mutate(start = pos - k, end = pos + k) %>%
    mutate(region = sprintf("%s:%i-%i", chr, start, end)) %>%    
    select(region) %>%
    write.table(., file = tmp.filename, sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

checkLowerCase <- function(char) {
  char == tolower(char)
}