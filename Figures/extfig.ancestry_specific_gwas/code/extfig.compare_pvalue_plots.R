# compare genome-wide significant signal (p<5e-8) from 12 meta to noAFR and to noASN
# GWAMA results without AFR and ASN cohorts

library(data.table)
library(ggplot2)

############

#meta12 vs noAFR

meta12_noAFR_compare <- fread("data_/AD_GWAMA_meta12_noAFR_compare_lead_snps.txt")
head(meta12_noAFR_compare)

library(ggplot2)
#p-value
ggplot(meta12_noAFR_compare, aes(x=meta12_noAFR_compare$`noAFR_-log10pvalue`, y=meta12_noAFR_compare$`12meta_-log10pvalue`)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  geom_abline(intercept=0, slope=1) +
  ggtitle("12 cohort meta-analysis p-value: with and without AFR cohorts")

ggsave(filename = "extfig.noAFR_v_meta12_pvalue.pdf", device = "pdf", width = 5, height = 5, units = "in")

#meta12 vs noASN
meta12_noASN_compare <- fread("data_/AD_GWAMA_meta12_noASN_compare_lead_snps.txt")
head(meta12_noASN_compare)

#p-value
ggplot(meta12_noASN_compare, aes(x=meta12_noASN_compare$`noASN_-log10pvalue`, y=meta12_noASN_compare$`12meta_-log10pvalue`)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  geom_abline(intercept=0, slope=1) +
  ggtitle("12 cohort meta-analysis p-value: with and without ASN cohorts")

ggsave(filename = "extfig.noASN_v_meta12_pvalue.pdf", device = "pdf", width = 5, height = 5, units = "in")

