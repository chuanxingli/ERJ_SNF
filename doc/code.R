library(DiagrammeR)


grViz("
digraph boxes_and_circles {
      
      # a 'graph' statement
      graph [overlap = true,style=filled, fontsize = 12,fontname = Helvetica]
      label = 'workflow'
      labelloc = 't'
      fontsize = 28
      # fontsize of label

      
      # several 'node' statements
      node [shape = box,
      fontsize = 12,
      fontname = Helvetica
      fillcolor = blue
      style = filled]
      A;B; C; D; E; F
      
      node [shape = circle,
      fixedsize = true,
      width = 0.9
      fillcolor=white] // sets as circles
      1; 2; 3; 4; 5; 6; 7; 8000
      
      node [shape=box]
      djf [label='ddd']

      node [shape=diamond]
      ste


      node [shape=parallelogram,fixedsize=F]
      IO [label='Input/Output']

     node [shape=record, width=1];
      b  [label='Process\\none\\n two three\\nfour five six seven\\n']
      d  [label='one\\l two three\\lfour five six seven\\l']
      f  [label='one\\r two three\\rfour five six seven\\r']

     node [shape = cylinder]
      g [label = 'storage']

     node [shape = note]
      h [label = 'note']

     node [shape = hexagon,
     color=red]
     i [label = 'preparation/\\nmanual work']

     node [shape = plaintext]
     j [label = 'start/\\end']
     

      # several 'edge' statements
      A->1 B->2 B->3 B->4 C->A
      1->D E->A 2->4 1->5 1->F
      E->6 4->6 5->7 6->7 3->8000 1->djf
      A->{4 2}
      b -> d -> f-> g-> h-> i-> j


      subgraph cluster {
           label = 'subgroup';
           fontsize = 14;
           rank = same; D; F; djf;
           D -> F -> djf [color = grey arrowhead=none]}
}",width=500,height=800)


# http://www.di.fc.ul.pt/~jpn/r/GraphicalTools/diagrammeR.html
# 
# download.file("https://wcs.smartdraw.com/flowchart/img/basic-symbols.jpg?bn=1510011143",destfile="~/Documents/R/flowchart.jpg")
# 
# # more details in https://www.smartdraw.com/flowchart/flowchart-symbols.htm
# 
# Node shape for graphviz
# https://www.graphviz.org/doc/info/shapes.html