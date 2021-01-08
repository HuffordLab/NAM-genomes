library(ggplot2)
color_fill <- c("red","darkred", "blue", "darkblue", "yellow", "gold", "green", "darkgreen")
a1 <- col2rgb(color_fill)
a2 <- rgb2hsv(a1)
n <- 0.4
color_x_axis <- c("gold", "blue", "blue", "blue", "blue", "blue", "blue", "grey", "grey", "grey", "pink", "orange", "orange", "green", "green", "green", "green", "green", "green", "green", "green", "green", "green", "green", "green", "green")
new_names <- c("B73", "B97", "Ky21", "M162W", "Ms71", "Oh43", "Oh7B", "M37W", "Mo18W", "Tx303", "HP301", "P39", "II14H", "CML52", "CML69", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "Ki3", "Ki11", "NC350", "NC358", "Tzi8")
P <- ggplot(all_NAM_file, aes(fill = factor(method, levels=c("A", "B", "C", "D", "E", "F", "G", "H")), y=Value, x=NAM)) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values = hsv(a2[1,], a2[2,]*n, a2[3,]),  breaks=c("A", "B", "C", "D", "E", "F" , "G" , "H"),labels=c("Ab initio Viridiplanteae", "Evidence Viridiplanteae", "Ab initio Poaceae", "Evidence Poaceae", "Ab initio Andropogoneae", "Evidence Andropogoneae", "Ab initio Maize", "Evidence Maize") ) + ylab("Number of gene models") + scale_x_discrete(labels= new_names) + scale_y_continuous(expand = c(0,0)) 
final_plot <- P + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = color_x_axis, size=20), axis.title=element_text(size=24), axis.text.y = element_text(size = 20), legend.text=element_text(size=20))
png("phylo.png", width=14,
    height=10,
    units="in",
    res=300)
print(final_plot)
dev.off()
