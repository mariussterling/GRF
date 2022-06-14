
library('glue')
df = read.table(glue("familysize_uniform",grids,".out"))
#df = read.table(glue("familysize_",grids,".out"))
load(file = glue("critical_val_",grids,".Rdata"))



df$INCOMED = df$INCOMED / 1000
df$ub = df$TAU + sqrt(df$VAR) * qnorm(0.975)
df$lb = df$TAU - sqrt(df$VAR) * qnorm(0.975)

df$uub_18 = df$TAU + sqrt(df$VAR) * uniform_q_star_18
df$ulb_18 = df$TAU - sqrt(df$VAR) * uniform_q_star_18

df$uub_22 = df$TAU + sqrt(df$VAR) * uniform_q_star_22
df$ulb_22 = df$TAU - sqrt(df$VAR) * uniform_q_star_22

average_width_18 = mean((df$uub_18-df$ulb_18))
average_width_22 = mean((df$uub_18-df$ulb_22))

#png(glue("new_plot_figure3_",grids,".png"), width = 900, height = 400)
png(glue("plot_figure3_uniform_",grids,".png"), width = 900, height = 400)
par(mfrow = c(1, 2), bg = "transparent")
p.df = df[df$AGEFSTM == 18, ]
plot(p.df$INCOMED, p.df$TAU,
     # ylim = range((p.df$lb), (p.df$ub),
     #         (p.df$ulb_18), (p.df$uub_18)),
     ylim = c(-0.75,0.75), 
     col = "black",
     type = "l",
     ylab = "CATE",
     xlab = "Father's Income [$1k/year]",
     sub = " Mother 18 years old at first birth")
lines(p.df$INCOMED, p.df$lb, type = "b",  lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ub, type = "b", lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ulb_18, lty = 3, col = "darkblue" )
lines(p.df$INCOMED, p.df$uub_18, lty = 3, col = "darkblue")
abline(h=0, lty = 3)

p.df = df[df$AGEFSTM == 22, ]
plot(p.df$INCOMED, p.df$TAU,
     # ylim = range((p.df$lb), (p.df$ub),
     #             (p.df$ulb_22), (p.df$uub_22)),
     ylim = c(-0.5,0.5),  
     col = "black",
     type = "l",
     ylab = "CATE",
     xlab = "Father's Income [$1k/year]",
     sub = "Mother 22 years old at first birth")
lines(p.df$INCOMED, p.df$lb, type = "b", lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ub, type = "b", lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ulb_22, lty = 3, col = "darkblue")
lines(p.df$INCOMED, p.df$uub_22, lty = 3, col = "darkblue")
abline(h=0, lty = 3)
dev.off()




df$INCOMED = df$INCOMED / 1000
df$ub = df$TAU + sqrt(df$VAR) * qnorm(0.975)
df$lb = df$TAU - sqrt(df$VAR) * qnorm(0.975)

df$uub_18 = df$TAU + sqrt(df$VAR) * uniform_q_star_18
df$ulb_18 = df$TAU - sqrt(df$VAR) * uniform_q_star_18

df$uub_22 = df$TAU + sqrt(df$VAR) * uniform_q_star_22
df$ulb_22 = df$TAU - sqrt(df$VAR) * uniform_q_star_22

average_width_18 = mean((df$uub_18-df$ulb_18))
average_width_22 = mean((df$uub_18-df$ulb_22))

#png(glue("plot_figure3_",grids,".png"), width = 900, height = 400)
png(glue("plot_figure3_uniform_",grids,".png"), width = 900, height = 400)
par(mfrow = c(1, 2), bg = "transparent")
p.df = df[df$AGEFSTM == 18, ]
plot(p.df$INCOMED, p.df$TAU,
    # ylim = range((p.df$lb), (p.df$ub),
     #         (p.df$ulb_18), (p.df$uub_18)),
    ylim = c(-1.5,1.5), 
    col = "black",
     type = "l",
     ylab = "CATE",
     xlab = "Father's Income [$1k/year]",
     sub = " Mother 18 years old at first birth")
lines(p.df$INCOMED, p.df$lb, type = "b",  lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ub, type = "b", lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ulb_18, lty = 3, col = "darkblue" )
lines(p.df$INCOMED, p.df$uub_18, lty = 3, col = "darkblue")
abline(h=0, lty = 3)

p.df = df[df$AGEFSTM == 22, ]
plot(p.df$INCOMED, p.df$TAU,
    # ylim = range((p.df$lb), (p.df$ub),
     #             (p.df$ulb_22), (p.df$uub_22)),
    ylim = c(-1.5,1.5),  
    col = "black",
     type = "l",
     ylab = "CATE",
     xlab = "Father's Income [$1k/year]",
     sub = "Mother 22 years old at first birth")
lines(p.df$INCOMED, p.df$lb, type = "b", lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ub, type = "b", lty = 3, col = "deeppink", cex=0.8)
lines(p.df$INCOMED, p.df$ulb_22, lty = 3, col = "darkblue")
lines(p.df$INCOMED, p.df$uub_22, lty = 3, col = "darkblue")
abline(h=0, lty = 3)
dev.off()

