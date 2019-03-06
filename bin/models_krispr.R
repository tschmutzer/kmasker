load("data_krispr.RData")
options(warn=-1)


## incoming data

args = commandArgs(trailingOnly=TRUE)

incoming <- data.frame(score=as.numeric(args[1]), score.1=as.numeric(args[2]), score.2=as.numeric(args[3]), score.3=as.numeric(args[4]), GC_all=as.numeric(args[5]), GC_dist=as.numeric(args[6]), GC_prox=as.numeric(args[7]), complexity=as.numeric(args[8]), ent=as.numeric(args[9]), take_score=as.numeric(args[10]), cover=as.numeric(args[11]))

## Make models

model_GC_all <- loess(eff ~ GC_all, data = data_c1, span=0.7)

model_comp_ent <- loess(eff ~ complexity + ent, data = data_c1, span=0.9)

model_GC_prox_dist <- loess(eff ~ GC_dist + GC_prox, data = data_c1, span=0.9)

combined_data <- data.frame(model_GC_prox_dist=predict(model_GC_prox_dist), model_comp_ent=predict(model_comp_ent), model_GC_all=predict(model_GC_all), eff=data_c1$eff)

model_combined <- loess(eff ~ model_GC_prox_dist + model_comp_ent + model_GC_all, data = combined_data, span=0.9)



## Predictions

pred_eff <- predict(model_combined, data.frame(model_GC_prox_dist=predict(model_GC_prox_dist, incoming[c(6,7)]),
											   model_comp_ent=predict(model_comp_ent, incoming[c(8,9)]),
											   model_GC_all=predict(model_GC_all, incoming[5])))


pred_eff <- ifelse(pred_eff < 0 , 0, pred_eff)
pred_eff <- ifelse(pred_eff > 1, 1, pred_eff)

											   
# output
print(unname(pred_eff))
