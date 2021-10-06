# Graph Theory Analysis Pipeline
Summary:


Theoretical Background: 

Traumatic brain injury (TBI) can be studied as a disorder of “disconnection” (citation needed). It often results in long-lasting and diffuse damages to white matter tracts (Johnson et al., 2013), and its impact on the functional organisation of the brain is also reflected in functional connectivity (FC) measures obtained with resting-state fMRI (rsMRI) (Xiao et al., 2015). Behaviourally, the sequelae of TBI can be relatively idiosyncratic, with symptoms ranging from physical (e.g., low energy) to cognitive (such as difficulties in concentration and in memory) to affective (e.g., with mood disturbances), all of which detrimental to quality of life (Smith et al., 2003). Graph theory analysis (GTA) has been used to characterise the connectomics changes in both TBI and depression (Caeyenberghs et al., 2017; Gong & He, 2015; Sharp et al., 2014; Yun & Kim, 2021). Compared to “low-level” measures examining functional connectivity emanating from a single region or between pairs of regions, graph theory can characterise topological properties for networks defined with FC, further incorporating information about how the “neighbours” of a region are connected, and potentially reveal functional integration and segregation across the network (Rubinov & Sporns, 2010). 

  

While ample work has been produced with GTA on FC data, the results have been relatively mixed for both TBI and for depression. Analysis of functional imaging data for TBI patients have sought to identify measures that distinguish patients from healthy controls, or measures that are associated with functional outcomes within the patient population. In general, TBI has been linked to alterations in overall FC level, as well as a more lattice-like organisation (Caeyenberghs et al., 2017). Specific to mood-related symptoms, one exploratory analysis, with GTA measures selected from the TBI literature, has found associations between depressive symptom and local graph theory metrics for nodes within DMN and FPCN (van der Horn et al., 2017). This observation corroborates the consistent link between primary depression and aberrant FC for midline versus lateral structures (Sundermann et al., 2014), as well as how TBI influences the DMN and induces changes in cognition and in executive functions (citation from behavioural bunch to be added). Motivated by the need for replication studies in neuroimaging (Fletcher & Grafton, 2013), the current investigation will test a depression-related hypothesis generated in van der Horn et al. (2017), namely that higher depressive scores are associated with higher measures of centrality in DMN-related regions and a reduction of centrality measures in executive network-related regions. 

  

Of note, the van der Horn et al. (2017) study focused on a cohort of mild TBI (mTBI) patients. On one hand, this prompts the question of whether said results could be generalised to TBI of greater severity. On the other hand, it has been reported that, compared to their family members, mTBI patients have more frequent health complaints, but the reverse is true for severe TBI patients (Seel et al., 2003), which may lead to a possible “inverted-U shape” in terms of self-reported depressive symptoms. It would hence be interesting to include general subjective complaints (e.g., as measured in SF36) as a covariate in further analysis, which may provide putative evidence for how “insight” – or symptom awareness – influences depression scores in TBI. 

Global properties 

Based on the premise that certain brain-wide network properties associated with major depressive disorder (MDD) may also underlie depressive symptoms co-morbid with other conditions, we would further investigate whether global graph theory metrics generally implicated in case-control studies for MDD could help predict the presence or absence of post-TBI depression. We would also consider measures that are likely with good interpretability (unlike measures based on path length), but less studied in the existing literature on MDD. 

  

As seen in recent reviews, there are mixed findings from GTA analysis on rsMRI for primary depression patients, which has been suggested to be attributable to the heterogeneity of the clinical samples (Gong & He, 2015; Yun & Kim, 2021). In terms of functional segregation, a variety of alterations in community structure have been reported (Yun & Kim, 2021), where a “community” describes a module of regions that are more tightly connected with each other than those outside the module. Globally, modularity score may increase for medication-naïve MDD patients (Ye et al., 2015), and novel brain stimulation or pharmacological interventions result in lower modularity in MDD patients (Caeyenberghs et al., 2018; Daws et al., 2021) – noting that Caeyenberghs et al. found no association between HDRS change and modularity change, but Daws did with BDI. Meanwhile, although TBI patients can also display disturbances to modular organisation, there seems to be no consensus on the directionality of effect on global modularity (Caeyenberghs et al., 2017), with one study reporting negative correlation between modularity and severity of post-concussion symptoms (measured with Rivermead Post-concussion Symptoms Questionnaire) six months after injury (Messé et al., 2013), though there may also be specific reductions in inter-modular connectivity (Han et al., 2014, 2016). Hence, research on TBI also warrants further investigation into how modular organisation contributes to the functional outcomes, and it may be necessary to examine participation coefficients (PC). This second measure would help characterise whether nodes exhibit more intra-modular (provincial hubs) or inter-modular (connector hubs) connections, and, when calculated for each node, has been linked to whether lesion to each region would produce more specific or broader behavioural deficits (Warren et al., 2014). 

  

In measuring functional integration, path length-based metrics such as global efficiency have potential interpretability issues regarding its implausible assumption on information propagation, though they do demonstrate relevance to behavioural variables (Fornito et al., 2013). Indeed, reduction in global efficiency was relatively characteristic of TBI (citation needed), and training on cognitive strategies has been reported to increase this metric whilst improving trail-making performance (Han et al., 2020). However, research on MDD reported contradictory findings (Meng et al., 2014; Zhang et al., 2011). Yun and Kim (2021) suggested that some of the contradictions may be due to the age of onset (e.g., adolescent or late-life) or illness durations that different studies focused on. Since emotional processing can undergo changes across age on both behavioural and neural levels (e.g., Schweizer et al., 2019), it is possible that the FC associated with depressive symptoms would also differ with age. Though the effect of age on functional connectivity for depression has not been systematically tested in MDD, we nevertheless plan to conduct an additional analysis to include age as a covariate or even a moderating factor, which might provide preliminary evidence for how age contributes to the heterogeneity with which neural changes can subserve depressive symptoms. 

  

Beyond segregation and integration, the topological properties of brain networks can be further evaluated with assortativity, which denotes the extent to which nodes preferentially connect to nodes of a similar degree. Functional connectivity networks seem to display positive assortativity (Eguíluz et al., 2005), and, since high-degree nodes (hubs) link to each other, a “rich-club” organisation would emerge, though this latter feature has more commonly been studied in the context of structural networks (Heuvel & Sporns, 2011). Importantly, given the theoretical links between assortativity and network resilience (Rubinov & Sporns, 2010), it may be able to capture TBI-related alterations in functional integrity beyond modularity and global efficiency, where the brain-wide but nevertheless selective damages to network hubs could be captured by lowered assortativity. To our knowledge, comparatively much fewer rsMRI-based studies have examined this third branch of network topology in TBI or in depression, though we are aware that this may be due to publication biases. Though one study has previously reported a null result between assortativity and post-TBI complaints not specific to affective symptoms (Messé et al., 2013), non-TBI studies have reported tentative links between lowered assortativity and primary depression or, at symptom level, suicide ideation (Chen et al., 2021; Erguzel et al., 2019; Lin et al., 2020; Wagner et al., 2019), together with lowered rich-club levels (Zhang et al., 2018). We would therefore include assortativity as a global metric in our analysis. 

 

Network-specific properties 

Functional connectivity changes under TBI have often been summarised and conceptualised in terms of large-scale functional networks, in particular the default mode network (DMN), the frontoparietal control network (FPCN, sometimes termed the central executive network, CEN), and the salience network (SN) (as reviewed in Hayes et al., 2016; Sharp et al., 2014; van der Horn et al., 2016). The same trio of networks have also been implicated in resting-state fMRI studies for major depressive disorder (MDD), alongside limbic regions (Kaiser et al., 2015; Mulders et al., 2015; Wang et al., 2012; Yun & Kim, 2021). An “imbalance” in network properties of DMN and FPCN nodes has been associated with depressive symptoms in mild TBI (van der Horn et al., 2017), but it is unknown whether the result would be generalisable to moderate and severe TBI patients. Meanwhile, in other neurological conditions such as strokes, severity depressive symptoms are positively correlated with connectivity between anterior DMN and the salience network (Balaev et al., 2018), though these alterations may not be consistently observed across studies (Liang et al., 2020). 

These results may be interpreted in light of the putative functional roles of DMN and FPCN (insert two possible reviews here), along with evidence for the salience network mediating their balance (Sridharan et al., 2008). Aberrant DMN activity in depression could be related to increase in internally-directed thinking and ruminations (Sheline et al., 2009; Zhou et al., 2020), while FPCN connectivity may be associated with variations in capacity for cognitive control (Seeley et al., 2007) and emotional regulation (Cole et al., 2014). However, the directionality of connectivity changes is not fully established for both TBI and primary depression, and previous studies do not always distinguish between inter- and intra-network FCs.  


Network Construction 

Nodes are defined by combining two data-driven parcellations of the cortex and the subcortical regions. The cortex is divided into 400 nodes (Schaefer et al., 2018), together with 54 subcortical nodes (Tian et al., 2020). For significant results, further analysis would seek to reproduce them in parcellations with coarser granularities (e.g., 100 cortical nodes and 16 subcortical nodes). It is known that edge definition would influence GTA outcomes for TBI (Han et al., 2016), where full Pearson’s correlation points to TBI-related hypoconnectivity despite partial correlations yielding the exact opposite, and more sophisticated approaches involving “graphical lasso” has been recommended for TBI research (Caeyenberghs et al., 2017). In view of future work addressing the serotonergic system (where using partial correlations that readily remove third-party influences would complicate an analysis designed to understand a third-party influence), we will use full correlations (Q: any corrections for autocorrelations). The outcome of this step would be a weighted but undirected graph with both positive and negative weights, which could be used to calculate graph theory metrics without thresholding (Rubinov & Sporns, 2011). Weight-conserving measures with negative weights will be especially useful for revealing the importance of anticorrelations between DMN and FPCN. 

  

Thresholding will be conducted to accommodate measures with sparsity requirements. Binary graphs will also be constructed at various thresholds, and the area under the metric-threshold curve (AUC) will be used to effectively pool the information from the range of thresholds covered. Weight distributions for depressed versus non-depressed patients (BDI >= 10 or < 10) will be plotted and compared visually prior to thresholding, and these low-level differences will be controlled for in further analysis. Taking a density-based thresholding (proportional thresholding) approach, the strongest 15%-35% connections (in steps of 1%) for each participant. will be included (though an algorithm-based method for determining the threshold will be used if available). 

 

For module decomposition, we applied the weight-conserving version of the Louvain algorithm (Rubinov & Sporns, 2011), which has previously been applied to patients of major depression (Lord et al., 2012). We used an asymmetric treatment for positive and negative weights. For fine-tuning, the modules obtained were used as the initial partition for another iteration of the algorithm, and this process continues until the modularity score reaches a maximum. The spatial resolution parameter, γ, has so far been set to default, though we are aware that there has been work trying to optimise its selection (Betzel et al., 2017). To mitigate the non-deterministic nature of modularity maximisation, and to merge the degenerate partitions with similar modularity scores into a single partition for each participant, we applied consensus clustering over the participant-specific partition ensemble obtained from 1000 repeats of the algorithm (Lancichinetti & Fortunato, 2012). This step requires a threshold, tau, over the agreement matrix, which was chosen to be 0.4 (Fornito et al., 2016), as the thresholds below 0.4 are generally well-performing for the Louvain algorithm, though with the caveat that, as it approaches zero, the final outcome would be increasingly divergent from the original decompositions present in the partition ensemble (Betzel et al., 2013). 

