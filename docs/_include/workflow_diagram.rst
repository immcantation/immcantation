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

       %% --- Entry point nodes (parallelogram = data/I-O symbol) ---
       raw[/"<i class='fa fa-file-text-o'></i>&nbsp;Raw Bulk FASTQ/FASTA"/]
       assembled[/"<i class='fa fa-file-text-o'></i>&nbsp;Assembled FASTA (Bulk, Single cell)"/]
       airr[/"<i class='fa fa-table'></i>&nbsp;AIRR-C TSV (10x)"/]

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

       subgraph tree["Lineage Tree Inference"]
           dowser["<b>Dowser</b>"]
       end       

       subgraph repertoire["Repertoire Analysis"]
           direction LR
           alakazam["<b>Alakazam</b>"]
           shazam["<b>SHazaM</b>"]
           %% invisible links keep Alakazam, SHazaM and Dowser side-by-side
           %% without drawing arrows between them
           alakazam ~~~ shazam
       end

       subgraph genot["Novel Allele, Genotyping&nbsp;"]
            %% &nbsp; padding widens the TIgGER node so the subgraph box is
            %% wide enough to contain the longer subgraph title without clipping
            tigger["<b>&nbsp;&nbsp;&nbsp;&nbsp;TIgGER&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</b>"]
       end
       subgraph embed["Embeddings"]
            amulety["<b>Amulety</b>"]
       end

       subgraph leg["Legend"]
           direction LR
           leg_data[/"Input data"/]
           leg_tool["Immcantation package"]
       end

       classDef stage fill:#dce9f5,stroke:#2980b9,stroke-width:1.5px,color:#2c3e50,rx:6,ry:6
       %% dashed border marks nodes as data files rather than tools
       classDef input fill:#ffffff,stroke:#888888,stroke-width:1.5px,color:#2c3e50,stroke-dasharray:4 2
       class presto,changeo,tigger,scoper,alakazam,shazam,dowser,amulety,leg_tool stage
       class raw,assembled,airr,leg_data input

       %% entry point connections
       raw --> pre  --> assembled
       assembled --> annot
       airr --> annot
       airr --> clonal

       annot --> clonal --> tree
       tree --> repertoire
       annot <--> genot --> clonal
       annot --> embed
       %% invisible link prevents Mermaid from routing an arrow between
       %% clonal and genot, which are not directly connected
       clonal ~~~ genot
