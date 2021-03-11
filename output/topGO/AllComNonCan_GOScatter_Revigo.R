# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0007043","cell-cell junction assembly",0.027,1.266,-7.007,3.790,-3.137,0.856,0.000),
c("GO:0097120","receptor localization to synapse",0.008,5.214,-2.374,3.269,-2.642,0.943,0.000),
c("GO:1905879","regulation of oogenesis",0.001,-3.076,2.808,2.314,-3.638,0.841,0.000),
c("GO:0007314","oocyte anterior/posterior axis specification",0.003,5.189,4.579,2.848,-2.975,0.728,0.004),
c("GO:0046710","GDP metabolic process",0.042,3.075,-0.396,3.979,-2.519,0.964,0.004),
c("GO:2000241","regulation of reproductive process",0.091,-5.515,1.427,4.313,-3.092,0.875,0.132),
c("GO:1904580","regulation of intracellular mRNA localization",0.001,-2.807,0.713,2.324,-1.821,0.889,0.132),
c("GO:0046716","muscle cell cellular homeostasis",0.007,-1.784,1.558,3.180,-1.761,0.892,0.137),
c("GO:0032486","Rap protein signal transduction",0.006,-5.655,-2.028,3.128,-2.818,0.855,0.147),
c("GO:0046956","positive phototaxis",0.000,-0.795,-3.441,1.778,-2.519,0.929,0.159),
c("GO:0048522","positive regulation of cellular process",1.903,-5.408,3.279,5.632,-2.951,0.803,0.224),
c("GO:1901700","response to oxygen-containing compound",0.537,-2.237,-3.088,5.082,-1.988,0.952,0.239),
c("GO:0016320","endoplasmic reticulum membrane fusion",0.002,1.665,-6.525,2.582,-2.642,0.883,0.252),
c("GO:0007028","cytoplasm organization",0.004,0.583,-6.832,2.948,-3.125,0.890,0.264),
c("GO:0097352","autophagosome maturation",0.008,0.007,-6.665,3.246,-1.821,0.884,0.274),
c("GO:0044087","regulation of cellular component biogenesis",0.447,-6.449,3.302,5.002,-1.895,0.855,0.279),
c("GO:0009966","regulation of signal transduction",0.849,-5.715,3.712,5.281,-2.336,0.748,0.299),
c("GO:0045175","basal protein localization",0.000,6.049,-2.040,1.756,-2.079,0.951,0.301),
c("GO:0051128","regulation of cellular component organization",0.988,-6.124,2.551,5.347,-2.253,0.846,0.304),
c("GO:0007169","transmembrane receptor protein tyrosine kinase signaling pathway",0.136,-5.905,-0.482,4.485,-2.338,0.814,0.305),
c("GO:0065008","regulation of biological quality",3.066,-5.898,2.277,5.839,-1.766,0.848,0.306),
c("GO:0007040","lysosome organization",0.021,1.457,-7.106,3.673,-1.917,0.857,0.312),
c("GO:0048518","positive regulation of biological process",2.070,-6.282,2.299,5.668,-2.796,0.844,0.313),
c("GO:0016336","establishment or maintenance of polarity of larval imaginal disc epithelium",0.000,5.797,3.421,1.380,-2.421,0.751,0.337),
c("GO:0043113","receptor clustering",0.008,5.874,-2.957,3.277,-2.342,0.909,0.357),
c("GO:0098586","cellular response to virus",0.010,-1.600,-3.259,3.339,-2.218,0.950,0.362),
c("GO:0030030","cell projection organization",0.670,1.147,-6.745,5.178,-2.465,0.852,0.367),
c("GO:1902875","regulation of embryonic pattern specification",0.000,-0.186,1.618,0.477,-1.761,0.894,0.372),
c("GO:0045167","asymmetric protein localization involved in cell fate determination",0.001,5.707,2.454,2.267,-1.866,0.772,0.412),
c("GO:0046940","nucleoside monophosphate phosphorylation",0.236,3.063,-0.493,4.725,-2.079,0.968,0.440),
c("GO:0000226","microtubule cytoskeleton organization",0.297,0.972,-7.091,4.825,-1.912,0.822,0.447),
c("GO:0120036","plasma membrane bounded cell projection organization",0.350,0.741,-6.530,4.896,-2.483,0.851,0.454),
c("GO:0009629","response to gravity",0.015,-2.740,-4.123,3.540,-1.780,0.944,0.473),
c("GO:0042332","gravitaxis",0.001,-1.262,-3.839,2.199,-1.780,0.926,0.495),
c("GO:0046598","positive regulation of viral entry into host cell",0.001,-1.877,4.561,2.516,-2.818,0.875,0.503),
c("GO:0045995","regulation of embryonic development",0.016,-5.165,0.823,3.552,-3.081,0.851,0.509),
c("GO:0043525","positive regulation of neuron apoptotic process",0.005,-3.707,5.378,3.029,-2.642,0.863,0.546),
c("GO:0007391","dorsal closure",0.003,5.436,4.545,2.780,-2.821,0.757,0.554),
c("GO:0000902","cell morphogenesis",0.311,5.173,4.919,4.845,-2.631,0.751,0.571),
c("GO:0016476","regulation of embryonic cell shape",0.001,-2.582,2.095,2.364,-1.800,0.864,0.580),
c("GO:0006930","substrate-dependent cell migration, cell extension",0.001,2.391,-5.566,2.371,-2.642,0.835,0.581),
c("GO:0007259","receptor signaling pathway via JAK-STAT",0.007,-6.020,-1.699,3.212,-1.725,0.842,0.587),
c("GO:0097696","receptor signaling pathway via STAT",0.007,-5.377,-1.700,3.214,-1.725,0.841,0.587),
c("GO:0007317","regulation of pole plasm oskar mRNA localization",0.001,-2.383,5.101,2.137,-1.821,0.857,0.592),
c("GO:0007279","pole cell formation",0.000,4.723,5.034,1.978,-2.042,0.782,0.595),
c("GO:0035162","embryonic hemopoiesis",0.012,5.158,4.331,3.422,-1.821,0.771,0.596),
c("GO:0060281","regulation of oocyte development",0.001,-2.992,3.105,2.243,-1.821,0.842,0.601),
c("GO:0032436","positive regulation of proteasomal ubiquitin-dependent protein catabolic process",0.027,-5.601,5.069,3.778,-1.708,0.849,0.623),
c("GO:0030714","anterior/posterior axis specification, follicular epithelium",0.000,5.635,3.796,0.602,-2.818,0.745,0.632),
c("GO:0035099","hemocyte migration",0.001,5.475,2.249,2.223,-1.917,0.744,0.640),
c("GO:0030710","regulation of border follicle cell delamination",0.004,-3.106,5.701,2.916,-3.119,0.840,0.641),
c("GO:0046627","negative regulation of insulin receptor signaling pathway",0.008,-5.256,4.193,3.271,-1.821,0.797,0.653),
c("GO:0048699","generation of neurons",0.320,4.816,4.567,4.858,-2.158,0.742,0.664),
c("GO:0000132","establishment of mitotic spindle orientation",0.014,4.015,-4.319,3.485,-1.891,0.729,0.670),
c("GO:0016332","establishment or maintenance of polarity of embryonic epithelium",0.000,5.771,3.626,1.643,-2.218,0.743,0.670),
c("GO:0046425","regulation of receptor signaling pathway via JAK-STAT",0.012,-5.342,4.430,3.442,-1.780,0.806,0.671),
c("GO:0051653","spindle localization",0.020,4.935,-3.498,3.643,-1.708,0.900,0.673),
c("GO:1904892","regulation of receptor signaling pathway via STAT",0.013,-5.367,3.329,3.473,-1.780,0.805,0.675),
c("GO:0001738","morphogenesis of a polarized epithelium",0.014,5.281,4.877,3.493,-2.688,0.769,0.681),
c("GO:0030713","ovarian follicle cell stalk formation",0.000,5.276,4.632,1.959,-2.007,0.764,0.687),
c("GO:0045196","establishment or maintenance of neuroblast polarity",0.001,5.475,3.415,2.220,-1.866,0.755,0.692));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

