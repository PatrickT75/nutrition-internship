library(rcompanion)

### one-way ANOVA of Shannon diversity vs kombucha, Tukey post-hoc
mydata = read.csv('bacterial alpha diversity.csv')
mydata$Kombucha = factor(mydata$Kombucha)
av = aov(mydata$Shannon ~ mydata$Kombucha)
tukey = TukeyHSD(av)
write.csv(tukey$`mydata$Kombucha`, 'alpha_bacterial_tukey.csv')

### show as compact-letter display
file = read.csv("alpha_bacterial_tukey.csv")
P.adj = file$p.adj
Comparison = file$comparison
cld = cldList(comparison = Comparison, p.value = P.adj)