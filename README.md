# kf2vec_data 
Data and scripts used for validation of kf2vec software at https://github.com/noraracht/kf2vec.

Populating...

## Raw data

* Models


  * Wol19 (with queries removed)
    - [Cladded/unchunked](https://github.com/noraracht/kf2vec_Wol19_Models/tree/main/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24)
    - [Uncladded/unchunked](https://github.com/noraracht/kf2vec_Wol19_Models/tree/main/k7_v37_8k_s28_TrainClassf_10K_TOL_Global)
    - [Cladded/chunked](https://github.com/noraracht/kf2vec_Wol19_Models/tree/main/k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks)
    - [Uncladded/chunked](https://github.com/noraracht/kf2vec_Wol19_Models/tree/main/k7_v37_8k_s28_TrainClassf_10K_TOL_Global_Chunks)
  
    For uncladded models, please `gunzip embeddings_subtree_0.csv.gz` before usage.

  * Wol19 (full)
    - [Cladded/unchunked]()
    - [Uncladded/unchunked]()


  * Wol23
    - [Cladded/chunked](https://github.com/noraracht/kf2vec_Wol23_Models/tree/main/k7_v42_8k_s28_train_model_65K_TOL_ChunksExp_MODEL)


  * Eukaryotes (cladded/unchunked)
    - [Fungi](https://github.com/noraracht/kf2vec_data/tree/main/k7_v62_8k_train_model_fungi_Claded_Unchunked)
    - [Insects](https://github.com/noraracht/kf2vec_data/tree/main/k7_v62_8k_train_model_insects_Claded_Unchunked_try2)
      

* Inputs k-mer frequencies and genomes (uplodad in progress)


## Results

<!---This section contains summary data tables and scripts we used to process them.--->

* Placement of full genomes
  
  * Model parameters, effects of divide-and-conquer, and chunking
    - The directory [full_genomes](https://github.com/noraracht/kf2vec_data/tree/main/full_genomes) contains results for the placement of full genomes. These summary tables serve as an input into the script [incomplete_genomes_git_full_genomes.R](https://github.com/noraracht/kf2vec_data/blob/main/full_genomes/incomplete_genomes_git_full_genomes.R) and were used to generate plots in Figures 2A, S1, S3, and S7.
      
  * Variable k-mer length
    - The directory [kmer_len](https://github.com/noraracht/kf2vec_data/tree/main/kmer_len) contains results for the variable k-mer length experiment, and the script [incomplete_genomes_git_kmerlen.R](https://github.com/noraracht/kf2vec_data/blob/main/kmer_len/incomplete_genomes_git_kmerlen.R) that was used to create Figure 2B.
  
  * Training progression
    - The directory [train_test_err](https://github.com/noraracht/kf2vec_data/tree/main/train_test_err) contains script [incomplete_genomes_train_test_err.R](https://github.com/noraracht/kf2vec_data/blob/main/train_test_err/incomplete_genomes_train_test_err.R) and data to generate diagrams in Figures S2A and S2B.


* Benchmarking
   
  * Comparison to CAFE
    - The directory [cafe_dist_comparison](https://github.com/noraracht/kf2vec_data/tree/main/cafe_dist_comparison) contains the summary results along with the script [placement_err_cafe_vs_us_git.R](https://github.com/noraracht/kf2vec_data/blob/main/cafe_dist_comparison/placement_err_cafe_vs_us_git.R) used to generate plots shown in Figures 2C and 2D.

  * Comparison to DEPP
    - The repository [depp](https://github.com/noraracht/kf2vec_data/tree/main/depp) contains data and the script [placement_err_qual_metrics_git_depp.R](https://github.com/noraracht/kf2vec_data/blob/main/depp/placement_err_qual_metrics_git_depp.R) to produce diagrams in Figures 2E, 2F, S4, S5 and S10.

  * Comparison to EPA-ng
    - Repository contains [fungi_results](https://github.com/noraracht/kf2vec_data/tree/main/k7_v62_8k_train_model_fungi_Claded_Unchunked_CmpClade) results on a fungal dataset that was used for comparison.
        
  * Comparison to Skmer
    - Repository contains [insects_results](https://github.com/noraracht/kf2vec_data/tree/main/k7_v62_8k_train_model_insects_Claded_Unchunked_CmpClade_try2) results on a insect dataset that was used for evaluation.


* Placement of incomplete genomes
  
  * Controlled length fragments
    - The [controlled_fragments](https://github.com/noraracht/kf2vec_data/tree/main/controlled_fragments) directory contains multiple summary result files that serve as inputs to the script [incomplete_genomes_git_fragments.R](https://github.com/noraracht/kf2vec_data/blob/main/controlled_fragments/incomplete_genomes_git_fragments.R), which was used to produce the diagrams in Figures 3A and S6.
  
  * Contigs
    - Summary results for contig placement are provided in the file [summary_TOL_full_genome_contigs_six_models.pl_error.gz](https://github.com/noraracht/kf2vec_data/blob/main/tol_contigs/summary_TOL_full_genome_contigs_six_models.pl_error.gz), which should be unarchived using `gunzip`. The resulting table serves as input to the script [placement_err_dev_queries_git_contigs.R](https://github.com/noraracht/kf2vec_data/blob/main/tol_contigs/placement_err_dev_queries_git_contigs.R), which generates the diagrams shown in Figures 3B and S8–9.
    
  * HiFi reads
    - Summary results for long-read placement are provided in the files [summary_hifi_six_models.pl_error_part1.csv.gz](https://github.com/noraracht/kf2vec_data/blob/main/hifi/summary_hifi_six_models.pl_error_part1.csv.gz) and [summary_hifi_six_models.pl_error_part2.csv.gz](https://github.com/noraracht/kf2vec_data/blob/main/hifi/summary_hifi_six_models.pl_error_part2.csv.gz), which should be unarchived using `gunzip`. hese files serve as input to the script [Summary_reads_placement_git_hifi.R](https://github.com/noraracht/kf2vec_data/blob/main/hifi/Summary_reads_placement_git_hifi.R), which produced Figures 3C and 3D.
    

    
* CAMI2
  - The directory [cami2](https://github.com/noraracht/kf2vec_data/tree/main/cami2) contains the CAMI2 analysis results generated by [AMBER](https://github.com/CAMI-challenge/AMBER). The script [placement_err_dev_queries_git_cami2.R](https://github.com/noraracht/kf2vec_data/blob/main/cami2/placement_err_dev_queries_git_cami2.R) was used to generate the plots shown in Figures 4A and 4B.
  
