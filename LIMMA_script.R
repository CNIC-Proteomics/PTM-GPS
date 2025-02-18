library(limma)
setwd('S:/U_Proteomica/LABS/LAB_BI/Proyecto_IR_PTMs/Ischemic_Analysis_JMR_withNextflow/ModStats-Navigator/tmp-vct227p2')

integrations <- c('Z_pgm2p', 'Z_pgm2p_dNM', 'Z_p2qf', 'Z_qf2q', 'Z_q2all')

intgr <- integrations[1]

groups <- c('20min', '40min', '80min', '120min', '6h', '12h', '24h', '120NOR', '24PreC')

for (g in groups) {
  for (intgr in integrations) {
    df <- read.csv(paste0('S:/U_Proteomica/LABS/LAB_BI/Proyecto_IR_PTMs/Ischemic_Analysis_JMR_withNextflow/ModStats-Navigator/tmp-vct227p2', '/', g, '/', intgr, '/', 'input.tsv'), header=T, sep='	')
    rownames(df) <- df$LEVEL
    df$LEVEL <- NULL

    samples <- colnames(df)
    design <- matrix(rep(1, length(samples)))
    rownames(design) <- samples
    colnames(design) <- g

    fit <- lmFit(df,design)
    fit <- eBayes(fit)

    pvalue <- as.data.frame(fit$p.value)
    meanRow <- as.data.frame(fit$Amean)

    output <- merge(meanRow, pvalue, 'row.names', all=T)
    colnames(output) <- c('LEVEL', 'Mean', 'pvalue')

    write.table(output, paste0('S:/U_Proteomica/LABS/LAB_BI/Proyecto_IR_PTMs/Ischemic_Analysis_JMR_withNextflow/ModStats-Navigator/tmp-vct227p2', '/', g, '/', intgr, '/', 'output.tsv'), sep='	', row.names=F, col.names=T)
  }
}
