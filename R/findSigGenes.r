

findSigGenes <- function(Expr, Label, Method = "tTest", Directed = TRUE, FdrCut = 0.01, FDCut = 1, SAMfdr = 0.20){
## for a given Expression profile, use statistic approach to evaluate the significance of difference between case and control samples;
## Expr: a data frame, the expression profile to be dealed with, the rownames should be the IDs of genes;
## Method: a string, should be one of the "tTest", "foldChange" or "SAM";
## Label: a vector of 0/1s, indicating the class of samples in the expression profile, 0 represents case,1 represents control;
## Directed: logical, if the the up or down regulated set should be distinguished;
## FdrCut: the fdr cutoff for T test;
## FDCut: the cutoff for fold change;
## SAMfdr: the parameter fdr.output used in SAM;
        if(Method == "tTest"){
                Tscore <- calT(Expr, Label);
                P_Val <- pt(abs(Tscore$T), lower.tail = FALSE, df = Tscore$df);
                Fdr <- p.adjust(P_Val);
                Result <- as.data.frame(cbind(rownames(Expr), Tscore$T, P_Val, Fdr));
                names(Result) <- c("ID", "t_statistic", "P", "FDR");
                save(Result, file = "tTestTable.rda");
                if(Directed == TRUE){
                        UpGeneSet <- Result[["ID"]][intersect(which(as.numeric(Result[["FDR"]]) < FdrCut), which(as.numeric(Result[["t_statistic"]]) > 0))];
                        LoGeneSet <- Result[["ID"]][intersect(which(as.numeric(Result[["FDR"]]) < FdrCut), which(as.numeric(Result[["t_statistic"]]) < 0))];
                        SigGenes <- list(UpGeneSet, LoGeneSet);
                        names(SigGenes) <- c("Up", "Low");
                        return(SigGenes);
                }else{
                        SigGenes <- Result[["ID"]][which(as.numeric(Result[["FDR"]]) < FdrCut)];
                        return(SigGenes);
                }
        }else if(Method == "foldChange"){
                Result <- calFD(Expr, Label);
                if(Directed == TRUE){
                        UpGeneSet <- Result[["ID"]][which(as.numeric(Result[["Log2FoldChange"]]) > FDCut)];
                        LoGeneSet <- Result[["ID"]][which(as.numeric(Result[["Log2FoldChange"]]) < (-FDCut))];
                        SigGenes <- list(UpGeneSet, LoGeneSet);
                        names(SigGenes) <- c("Up", "Low");
                        return(SigGenes);

                }else{
                        SigGenes <- Result[["ID"]][which(abs(as.numeric(Result[["Log2FoldChange"]])) > FDCut)];
                        return(SigGenes);
                }

        }else if(Method == "SAM"){
                        Label[which(Label == 0)] <- 2;
                        Result <- SAM(x = Expr, y = Label, resp.type = "Two class unpaired", geneid = rownames(Expr), fdr.output = SAMfdr);
                        UpGeneSet <- Result$siggenes.table$genes.up[,2];
                        LoGeneSet <- Result$siggenes.table$genes.lo[,2];
                        if(Directed == TRUE){
                                SigGenes <- list(UpGeneSet, LoGeneSet);
                                names(SigGenes) <- c("Up", "Low");
                                return(SigGenes);
                        }else{
                                SigGenes <- c(UpGeneSet, LoGeneSet);
                                return(SigGenes);
                        }
        }
}




