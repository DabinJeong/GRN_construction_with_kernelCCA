datasets = list(config['datasets'].keys())

results = []
for dataset in datasets:
        result_dat = []
        for i in range(1,int(config['datasets'][dataset]['num_tp'])):
                result_dat.append('result/{taxon}_{dataset}_{nComp}_norm{norm}_GRN_{threshold}_tp{i}.pkl'.format(taxon=config['datasets'][dataset]['taxon'],dataset=config['datasets'][dataset]['GEO_id'],nComp=config['datasets'][dataset]['num_comp'],norm=config['datasets'][dataset]['norm'],threshold=config['parameters']['edge_thr'],i=i))
        results.extend(result_dat)


rule all:
        input:
                results 

rule split_into_TF_nwk:
        input:
                TF_li = "TF_li/{taxon}_TF_list",
                PPI = "network_ppi/PPI_{taxon}.tsv"
        output:
                PPI_TF = "network_ppi/{taxon}_template_nwk_TF.tsv"
        run:
                import pandas as pd
                import numpy as np
                df_TF = pd.read_csv(input.TF_li,sep='\t')
                df_nwk = pd.read_csv(input.PPI,sep='\t') 
                df_nwk_TF = df_nwk.loc[lambda x:np.logical_and(x.protein1.isin(df_TF.Name),x.protein2.isin(df_TF.Name))]
                df_nwk_TF.to_csv(output.PPI_TF,sep='\t',index=False)               

rule split_into_TG_nwk:
        input:
                TF_li = "TF_li/{taxon}_TF_list",
                PPI = "network_ppi/PPI_{taxon}.tsv"
        output:
                PPI_TG = "network_ppi/{taxon}_template_nwk_TG.tsv"
        run:
                import pandas as pd
                import numpy as np
                df_TF = pd.read_csv(input.TF_li,sep='\t')
                df_nwk = pd.read_csv(input.PPI,sep='\t')
                df_nwk_TG = df_nwk.loc[lambda x:np.logical_and(~x.protein1.isin(df_TF.Name),~x.protein2.isin(df_TF.Name))]
                df_nwk_TG.to_csv(output.PPI_TG,sep='\t',index=False)

rule instantiate_TF_nwk:
        input:
                template_nwk_TF = "network_ppi/{taxon}_template_nwk_TF.tsv"
        params:
                corrCut = config["parameters"]["corr_thr"]
        output:
                inst_nwk_TF  = "TFTG_cluster/{taxon}_{dataset}_condition_specific_nwk_TF.tsv"
        shell:
                """
                python 0_instantiate_nwk.py {input.template_nwk_TF} data/{wildcards.taxon}_{wildcards.dataset}_gene_exp_profile  -corrCut {params.corrCut} -o {output.inst_nwk_TF} -nThreads 50
                """      

rule instantiate_TG_nwk:
        input:
                template_nwk_TG = "network_ppi/{taxon}_template_nwk_TG.tsv"
        params:
                corrCut = config["parameters"]["corr_thr"]        
        output:
                inst_nwk_TG  = "TFTG_cluster/{taxon}_{dataset}_condition_specific_nwk_TG.tsv"
        shell:
                """
                python 0_instantiate_nwk.py {input.template_nwk_TG} data/{wildcards.taxon}_{wildcards.dataset}_gene_exp_profile -corrCut {params.corrCut} -o {output.inst_nwk_TG} -nThreads 50
                """

rule cluster_TF_nwk:
        input: 
                inst_nwk_TF = "TFTG_cluster/{taxon}_{dataset}_condition_specific_nwk_TF.tsv"
        output: 
                cluster_TF = "TFTG_cluster/{taxon}_{dataset}_community_TF.tsv",
                network_TF = "TFTG_cluster/{taxon}_{dataset}_communityEdgelist_TF.tsv" 
        shell:
                """
                Rscript 1_clustering.R {input.inst_nwk_TF} {output.cluster_TF} {output.network_TF} 
                """

rule cluster_TG_nwk:
        input:
                inst_nwk_TG = "TFTG_cluster/{taxon}_{dataset}_condition_specific_nwk_TG.tsv"
        output:
                cluster_TG = "TFTG_cluster/{taxon}_{dataset}_community_TG.tsv",
                network_TG = "TFTG_cluster/{taxon}_{dataset}_communityEdgelist_TG.tsv"
        shell:
                """
                Rscript 1_clustering.R {input.inst_nwk_TG} {output.cluster_TG} {output.network_TG}
                """

rule cluster2dict_TF:
        input:
                cluster_TF = "TFTG_cluster/{taxon}_{dataset}_community_TF.tsv"
        output:
                dict_cluster_TF = "TFTG_cluster/{taxon}_{dataset}_community_TF.pkl"
        run:
                import pandas as pd
                import pickle
                df_TF_community_groupBy = pd.read_csv(input.cluster_TF,sep='\t')
                with open(output.dict_cluster_TF,'wb') as f:
                        pickle.dump(df_TF_community_groupBy.groupby('community')['gene'].apply(list).to_dict(),f)

rule cluster2dict_TG:
        input:
                cluster_TG = "TFTG_cluster/{taxon}_{dataset}_community_TG.tsv"
        output:
                dict_cluster_TG = "TFTG_cluster/{taxon}_{dataset}_community_TG.pkl"
        run:
                import pandas as pd
                import pickle
                df_TG_community_groupBy = pd.read_csv(input.cluster_TG,sep='\t')
                with open(output.dict_cluster_TG,'wb') as f:
                        pickle.dump(df_TG_community_groupBy.groupby('community')['gene'].apply(list).to_dict(),f)

rule inferGRN:
        input:
                dict_cluster_TF = "TFTG_cluster/{taxon}_{dataset}_community_TF.pkl",
                dict_cluster_TG = "TFTG_cluster/{taxon}_{dataset}_community_TG.pkl"
        params:
                nThreads = config["parameters"]["nThreads"],
                reg = config["parameters"]["reg_kernel"] 
        output:
                "result/{taxon}_{dataset}_{nComp}_norm{norm}_GRN_{thr}_tp{i}.pkl"
        shell:
                """
                python 2_GRN_inference.py network_GRN/GRN_{wildcards.taxon}.tsv data/{wildcards.taxon}_{wildcards.dataset}_gene_exp_profile.tp{wildcards.i} data/{wildcards.taxon}_{wildcards.dataset}_DEG.tp{wildcards.i} {output} -TFli TF_li/{wildcards.taxon}_TF_list -TFmodule {input.dict_cluster_TF} -TGmodule {input.dict_cluster_TG} -nComp {wildcards.nComp} -thr {wildcards.thr} -reg {params.reg} -nThreads {params.nThreads} -normalize {wildcards.norm}
                """


