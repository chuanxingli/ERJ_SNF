library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(svglite)
library(rsvg)
library(png)

## Date = 2019-03-07
semplot <- grViz("
digraph boxes_and_circles {

        # a 'graph' statement
        graph [overlap = true, style=filled, fontsize = 12,fontname = Helvetica, compound = true]


subgraph cluster_m {
        label = 'Design of analysis workflow'
        labelloc = 't'
        fontsize = 28
        # fontsize of label
        color = white
        fillcolor = white
                # several 'node' statements
                node []
                        x0 [label = 'Raw data',shape = cylinder]
                        x1 [label = 'Data processing\\n(process.R)',shape = box]

                        a [label='Processed data\\n(pdata.rdata)', shape = parallelogram, margin = 0]

                        

                subgraph cluster {
                        label = 'Unsupervised clustering\\n(SNFunsupervised.R)';
                        fillcolor = white;
                        color = black;
                        fontsize = 14;
                        g1 [label = 'Parameter optimization' shape = box]
                        g2 [label = 'Best cluster\\nnumbers estimation' shape = box] }

                subgraph cluster2 {
                        label = 'Supervised clustering\\n(SNFsupervised.R)';
                        fillcolor = white;
                        color = black;
                        fontsize = 14;
                        h1 [label = 'Permutation' shape = box]
                        h2 [label = 'Parameter optimization' shape = box]}


                node [shape = box]
                        d [label = 'Predicted clusters(pcluster.rdata)', shape = parallelogram, margin = 0]
                        e [label = 'Accuracy of prediction\\n(accuracy.R)',shape = box]
                        f [label = 'Accuracies\\n(acc.rdata)', shape = parallelogram, margin = 0]

                        i [label = 'Downstream Analysis']
                        # i1 [label = 'Blockcluster of \\nfunctions and communities']
                        # j [label = '3-nodes motifs']
                        # k [label = 'Candidate gene list \\nnetworks', shape = parallelogram, margin = 0]
                        # l [label = 'Classifier']
                        # m [label = 'Experimental Validation']
                        # n [label = 'Literature/DB \\nsupport' shape = hexagon]


                # several 'edge' statements
                x0 -> x1 -> a
                a -> g1 [lhead = cluster]
                a -> h1 [lhead = cluster2]
                g1 -> g2
                h1 -> h2
                g2 -> d [ltail = cluster]
                h2 -> d [ltail = cluster2]
                d -> e -> f-> i


}

        subgraph cluster_l {

                label = 'Legend';
                fontsize = 14;
                fillcolor = white;

                node [margin = 0.1, style = filled, fillcolor = grey]
                        ll1 [label = '', shape = cylinder]
                        ll2 [label = '', shape = box]
                        ll3 [label = '', shape = hexagon]
                        ll4 [label = '', shape = parallelogram]
                        ll5 [label = '', shape = diamond]

                node [shape = plaintext, style = unfilled]
                        ll1t [label = 'Data']
                        ll2t [label = 'Processing']
                        ll3t [label = 'Manually\\nwork',]
                        ll4t [label = 'Input/\\nOutput']
                        ll5t [label = 'Decision']
                #shape=parallelogram
                ll1 -> ll2 -> ll3 -> ll4 ->ll5[color=none style = invis]
                ll1 -> ll1t [color=none style = invis]
                ll2 -> ll2t [color=none style = invis]
                ll3 -> ll3t [color=none style = invis]
                ll4 -> ll4t [color=none style = invis]
                ll5 -> ll5t [color=none style = invis]
                {rank = same; ll1; ll1t}
                {rank = same; ll2; ll2t}
                {rank = same; ll3; ll3t}
                {rank = same; ll4; ll4t}
                {rank = same; ll5; ll5t}

}
        # align two subgraph
        newrank=true;
        {rank = same; ll5; i}


}",width=500,height=800)


filename <- paste0("design_",Sys.Date(),".pdf")
# semplot %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG(filename)
semplot %>% export_svg %>% charToRaw %>% rsvg_pdf(filename)
