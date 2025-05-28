#### v2.0.1

1. fix the bug of backgrounds for karyo plots 
2. fix the bug of link_colors for karyo text plot 

#### v2.0.4

1. fix the bug of tree plot 
   + when the branch length was less than 0, the dotted line will beyond the horizontal line
   + can't draw the tree when its mutil-way tree
2. add the arrow attribue for highlight block of circular karyo plot

#### v2.0.5

1. add the mow plot type for hcgd

#### v2.0.6

1. fix the bug of dividing function in STAT.pm

#### v2.0.7

1. add the 'chromosomes_stroke_width' options for ideogram of karyo plot 

#### v2.2.0

1. the unit of tree can be calc auto
2. add 'flower_plot' option for venn diagram (can draw flower plot with less 6 sets)
3. add two symbols (fan and semi-circle)

### v2.2.1

1. fix the unrooted tree layout bug (tree is away from the scale line)
    calc the tree size, then transform the tree group with translate and scale
    (can also use the 'getBBox' method in javascript, but is not support by inkscape)
2. add 'show_outrange' function for <ggplot2> type fig, default is 'no' (will clip the diagram by the xy coordinate)
3. add the gradient legend 

### v2.3.2

1. update the enrichPlot.pl script, add some functions 
2. add some options for LEGEND.pm and SYMBOL.pm, you can add the lines and text on the item
3. add the 'chromosomes_stroke_color' options for ideogram of karyo plot
