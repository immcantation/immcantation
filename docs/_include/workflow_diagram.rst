.. mermaid::

   %%{init: {
     'theme': 'base',
     'themeVariables': {
       'primaryColor': '#dce9f5',
       'primaryTextColor': '#2c3e50',
       'primaryBorderColor': '#2980b9',
       'lineColor': '#2980b9',
       'clusterBkg': '#f0f7fd',
       'clusterBorder': '#2980b9',
       'titleColor': '#2c3e50',
       'edgeLabelBackground': '#ffffff',
       'fontFamily': 'Lato, sans-serif',
       'fontSize': '16px'
     }
   }}%%
   flowchart LR

       subgraph pre["Pre-processing"]
           presto["<b>pRESTO</b>"]
       end

       subgraph annot["V(D)J Assignment"]
           direction TB
           changeo["<b>Change-O</b>"]
       end

       subgraph clonal["Clonal Clustering"]
           scoper["<b>SCOPer</b>"]
       end

       subgraph repertoire["Repertoire Analysis"]
           direction LR
           alakazam["<b>Alakazam</b>"]
           shazam["<b>SHazaM</b>"]
           dowser["<b>Dowser</b>"]
           %% invisible links keep Alakazam, SHazaM and Dowser side-by-side
           %% without drawing arrows between them
           alakazam ~~~ shazam
           alakazam ~~~ dowser
       end

       subgraph genot["Novel Allele, Genotyping&nbsp;"]
            %% &nbsp; padding widens the TIgGER node so the subgraph box is
            %% wide enough to contain the longer subgraph title without clipping
            tigger["<b>&nbsp;&nbsp;&nbsp;&nbsp;TIgGER&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</b>"]
       end
       subgraph embed["Embeddings"]
            amulety["<b>Amulety</b>"]
       end

       classDef stage fill:#dce9f5,stroke:#2980b9,stroke-width:1.5px,color:#2c3e50,rx:6,ry:6
       class presto,changeo,tigger,scoper,alakazam,shazam,dowser,amulety stage

       pre --> annot --> clonal --> repertoire
       annot <--> genot
       annot --> embed
       %% invisible link prevents Mermaid from routing an arrow between
       %% clonal and genot, which are not directly connected
       clonal ~~~ genot
