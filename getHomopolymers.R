getHomopolymers <- function(HRun.vcf, VCFPolyx.vcf) {
  Polyx.cmd <- sprintf("bcftools query -i 'POLYX>1' -f'%%CHROM:%%POS-%%INFO/POLYX\n' %s", VCFPolyx.vcf)
  Polyx.data <- system(Polyx.cmd, intern = TRUE)
  HRun.cmd <- sprintf("bcftools query -i 'HRun>1' -f'%%CHROM:%%POS-%%INFO/HRun\n' %s", HRun.vcf)
  HRun.data <- system(HRun.cmd, intern = TRUE)
  Polyx.data %>%
    as_tibble() %>%
    separate(value, into = c("id", "hp.length"), sep = "-") -> Polyx.data
  HRun.data %>%
    as_tibble() %>%
    separate(value, into = c("id", "hp.length"), sep = "-") -> HRun.data
  full_join(Polyx.data, HRun.data, by = "id") %>%
    rename(hp.length.Polyx = hp.length.x, hp.length.HRun = hp.length.y) %>%
    mutate(hp.length = pmax(hp.length.Polyx, hp.length.HRun, na.rm = TRUE),
           hp.length = as.numeric(hp.length),
           hp.source = if_else(is.na(hp.length.HRun), "Polyx", "HRun_or_Polyx")) %>%
    select(id, hp.length, hp.source)
}