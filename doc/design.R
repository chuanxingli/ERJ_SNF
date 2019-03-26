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
                        add0 [label = 'Parameters for random data\\nnomics, ngroup,\\nsample/gene sizes', shape = hexagon]
                        add1 [label = 'Generate random data\\n(generate_random_dataL.R)', shape = box]
                        x0 [label = 'Raw data',shape = cylinder]
                        x1 [label = 'Data processing\\n(generate_input.R)',shape = box]
                        a [label='Processed data\\n(input)', shape = parallelogram, margin = 0]
                       add2 [label='Similarity Network Fusion\\n(SNFtool)', shape = box]
                        add3 [label = 'Similarity Matrix\\n(Wi, W)', shape = parallelogram, margin = 0]


                subgraph cluster {
                        label = 'Unsupervised clustering\\n(analysis.R)';
                        fillcolor = white;
                        color = black;
                        fontsize = 14;
                        #g1 [label = 'Parameter optimization' shape = box]
                        g2 [label = 'Best cluster\\nnumbers estimation' shape = box]
                        g1 [label = 'Fixed\\ncluster number' shape = box]
                        g3 [label = 'Unsupervised clustering' shape = box]}

                subgraph cluster2 {
                        label = 'Supervised clustering\\n(analysis.R)';
                        fillcolor = white;
                        color = black;
                        fontsize = 14;
                        h1 [label = 'Permutation\\n(LOOCV)' shape = box]
                        h11 [label = 'Supervised prediction' shape = box]}

                subgraph cluster_3{
                        label = 'Generate Dataset';
                        fillcolor = white;
                        color = black;
                        fontsize = 14;
                        add0; add1
                }
                #{rank = same; add1; x0}

                node [shape = box]

                        d [label = 'Accuracy\\n(NMI)', shape = parallelogram, margin = 0]
                        # e [label = 'Parameter Optimization\\n(unupload.R)',shape = box]
                        # f [label = 'Accuracies\\n(resNMI.rdata)', shape = parallelogram, margin = 0]
                        # 
                        # i [label = 'Downstream Analysis']
                        # i1 [label = 'Blockcluster of \\nfunctions and communities']
                        # j [label = '3-nodes motifs']
                        # k [label = 'Candidate gene list \\nnetworks', shape = parallelogram, margin = 0]
                        # l [label = 'Classifier']
                        # m [label = 'Experimental Validation']
                        # n [label = 'Literature/DB \\nsupport' shape = hexagon]


                # several 'edge' statements
                x0 -> x1 -> a
                add3 -> g1 [lhead = cluster]
                add3 -> h1 [lhead = cluster2]
                {g1,g2} -> g3
                h1 -> h11
                g3 -> d [ltail = cluster]
                h11 -> d [ltail = cluster2]
                #d -> e -> f-> i
add0 -> add1 -> x0
a -> add2 [label = 'Parameters for SNF\\nK, alpha, t']
add2 -> add3


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
        {rank = same; ll5; d}


}",width=500,height=800)


filename <- paste0("./doc/design_",Sys.Date(),".png")
semplot %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG(filename)

