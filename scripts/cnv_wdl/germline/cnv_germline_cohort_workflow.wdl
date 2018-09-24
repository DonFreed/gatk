# Workflow for creating a GATK GermlineCNVCaller denoising model and generating calls given a list of normal samples. Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by PreprocessIntervals.padding (default 250)
#   and split into bins of length specified by PreprocessIntervals.bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the chromosomes of interest.
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_germline_cohort_workflow.wdl -i my_parameters.json
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCohortWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    File intervals
    File? blacklist_intervals
    Array[String]+ normal_bams
    Array[String]+ normal_bais
    String cohort_entity_id
    File contig_ploidy_priors
    Int num_intervals_per_scatter
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    # If true, AnnotateIntervals will be run to create GC annotations and explicit
    # GC correction will be performed by the model generated by
    Boolean? do_explicit_gc_correction
    File? gatk4_jar_override
    Int? preemptible_attempts

    ####################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    Int? padding
    Int? bin_length

    ##################################################
    #### optional arguments for AnnotateIntervals ####
    ##################################################
    File? mappability_track
    File? segmental_duplication_track
    Int? feature_query_lookahead
    Int? mem_gb_for_annotate_intervals

    ##############################################
    #### optional arguments for CollectCounts ####
    ##############################################
    String? collect_counts_format
    Int? mem_gb_for_collect_counts

    ########################################################################
    #### optional arguments for DetermineGermlineContigPloidyCohortMode ####
    ########################################################################
    Float? ploidy_mean_bias_standard_deviation
    Float? ploidy_mapping_error_rate
    Float? ploidy_global_psi_scale
    Float? ploidy_sample_psi_scale
    Int? mem_gb_for_determine_germline_contig_ploidy
    Int? cpu_for_determine_germline_contig_ploidy

    ############################################################
    #### optional arguments for GermlineCNVCallerCohortMode ####
    ############################################################
    Float? gcnv_p_alt
    Float? gcnv_p_active
    Float? gcnv_cnv_coherence_length
    Float? gcnv_class_coherence_length
    Int? gcnv_max_copy_number
    Int? mem_gb_for_germline_cnv_caller
    Int? cpu_for_germline_cnv_caller

    # optional arguments for germline CNV denoising model
    Int? gcnv_max_bias_factors
    Float? gcnv_mapping_error_rate
    Float? gcnv_interval_psi_scale
    Float? gcnv_sample_psi_scale
    Float? gcnv_depth_correction_tau
    Float? gcnv_log_mean_bias_standard_deviation
    Float? gcnv_init_ard_rel_unexplained_variance
    Int? gcnv_num_gc_bins
    Float? gcnv_gc_curve_standard_deviation
    String? gcnv_copy_number_posterior_expectation_mode
    Boolean? gcnv_enable_bias_factors
    Int? gcnv_active_class_padding_hybrid_mode

    # optional arguments for Hybrid ADVI
    Float? gcnv_learning_rate
    Float? gcnv_adamax_beta_1
    Float? gcnv_adamax_beta_2
    Int? gcnv_log_emission_samples_per_round
    Float? gcnv_log_emission_sampling_median_rel_error
    Int? gcnv_log_emission_sampling_rounds
    Int? gcnv_max_advi_iter_first_epoch
    Int? gcnv_max_advi_iter_subsequent_epochs
    Int? gcnv_min_training_epochs
    Int? gcnv_max_training_epochs
    Float? gcnv_initial_temperature
    Int? gcnv_num_thermal_advi_iters
    Int? gcnv_convergence_snr_averaging_window
    Float? gcnv_convergence_snr_trigger_threshold
    Int? gcnv_convergence_snr_countdown_window
    Int? gcnv_max_calling_iters
    Float? gcnv_caller_update_convergence_threshold
    Float? gcnv_caller_internal_admixing_rate
    Float? gcnv_caller_external_admixing_rate
    Boolean? gcnv_disable_annealing

    ###################################################
    #### arguments for PostprocessGermlineCNVCalls ####
    ###################################################
    Int ref_copy_number_autosomal_contigs
    Array[String]? allosomal_contigs

    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    if (select_first([do_explicit_gc_correction, false])) {
        call CNVTasks.AnnotateIntervals {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                mappability_track = mappability_track,
                segmental_duplication_track = segmental_duplication_track,
                feature_query_lookahead = feature_query_lookahead,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_annotate_intervals,
                preemptible_attempts = preemptible_attempts
        }
    }

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call CNVTasks.CollectCounts {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts
        }
    }

    call DetermineGermlineContigPloidyCohortMode {
        input:
            cohort_entity_id = cohort_entity_id,
            read_count_files = CollectCounts.counts,
            contig_ploidy_priors = contig_ploidy_priors,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_determine_germline_contig_ploidy,
            cpu = cpu_for_determine_germline_contig_ploidy,
            mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
            mapping_error_rate = ploidy_mapping_error_rate,
            global_psi_scale = ploidy_global_psi_scale,
            sample_psi_scale = ploidy_sample_psi_scale,
            preemptible_attempts = preemptible_attempts
    }

    call CNVTasks.ScatterIntervals {
        input:
            interval_list = PreprocessIntervals.preprocessed_intervals,
            num_intervals_per_scatter = num_intervals_per_scatter,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    scatter (scatter_index in range(length(ScatterIntervals.scattered_interval_lists))) {
        call GermlineCNVCallerCohortMode {
            input:
                scatter_index = scatter_index,
                cohort_entity_id = cohort_entity_id,
                read_count_files = CollectCounts.counts,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                intervals = ScatterIntervals.scattered_interval_lists[scatter_index],
                annotated_intervals = AnnotateIntervals.annotated_intervals,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_germline_cnv_caller,
                cpu = cpu_for_germline_cnv_caller,
                p_alt = gcnv_p_alt,
                p_active = gcnv_p_active,
                cnv_coherence_length = gcnv_cnv_coherence_length,
                class_coherence_length = gcnv_class_coherence_length,
                max_copy_number = gcnv_max_copy_number,
                max_bias_factors = gcnv_max_bias_factors,
                mapping_error_rate = gcnv_mapping_error_rate,
                interval_psi_scale = gcnv_interval_psi_scale,
                sample_psi_scale = gcnv_sample_psi_scale,
                depth_correction_tau = gcnv_depth_correction_tau,
                log_mean_bias_standard_deviation = gcnv_log_mean_bias_standard_deviation,
                init_ard_rel_unexplained_variance = gcnv_init_ard_rel_unexplained_variance,
                num_gc_bins = gcnv_num_gc_bins,
                gc_curve_standard_deviation = gcnv_gc_curve_standard_deviation,
                copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                enable_bias_factors = gcnv_enable_bias_factors,
                active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                learning_rate = gcnv_learning_rate,
                adamax_beta_1 = gcnv_adamax_beta_1,
                adamax_beta_2 = gcnv_adamax_beta_2,
                log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                min_training_epochs = gcnv_min_training_epochs,
                max_training_epochs = gcnv_max_training_epochs,
                initial_temperature = gcnv_initial_temperature,
                num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                max_calling_iters = gcnv_max_calling_iters,
                caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                disable_annealing = gcnv_disable_annealing,
                preemptible_attempts = preemptible_attempts
        }
    }

    Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCohortMode.gcnv_call_tars)

    scatter (sample_index in range(length(CollectCounts.entity_id))) {
        call CNVTasks.PostprocessGermlineCNVCalls {
            input:
                calling_configs = GermlineCNVCallerCohortMode.calling_config_json,
                denoising_configs = GermlineCNVCallerCohortMode.denoising_config_json,
                gcnvkernel_version = GermlineCNVCallerCohortMode.gcnvkernel_version_json,
                sharded_interval_lists = GermlineCNVCallerCohortMode.sharded_interval_list,
                entity_id = CollectCounts.entity_id[sample_index],
                gcnv_calls_tars = call_tars_sample_by_shard[sample_index],
                gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                allosomal_contigs = allosomal_contigs,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                sample_index = sample_index,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts
        }
    }

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals
        Array[File] read_counts_entity_ids = CollectCounts.entity_id
        Array[File] read_counts = CollectCounts.counts
        File contig_ploidy_model_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_model_tar
        File contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar
        Array[File] gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar
        Array[Array[File]] gcnv_calls_tars = GermlineCNVCallerCohortMode.gcnv_call_tars
        Array[File] gcnv_tracking_tars = GermlineCNVCallerCohortMode.gcnv_tracking_tar
        Array[File] genotyped_intervals_vcfs = PostprocessGermlineCNVCalls.genotyped_intervals_vcf
        Array[File] genotyped_segments_vcfs = PostprocessGermlineCNVCalls.genotyped_segments_vcf
    }
}

task DetermineGermlineContigPloidyCohortMode {
    String cohort_entity_id
    Array[File] read_count_files
    File contig_ploidy_priors
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    # Model parameters
    Float? mean_bias_standard_deviation
    Float? mapping_error_rate
    Float? global_psi_scale
    Float? sample_psi_scale

    # We do not expose Hybrid ADVI parameters -- the default values are decent

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=${default=8 cpu}
        export OMP_NUM_THREADS=${default=8 cpu}

        gatk --java-options "-Xmx${command_mem_mb}m"  DetermineGermlineContigPloidy \
            --input ${sep=" --input " read_count_files} \
            --contig-ploidy-priors ${contig_ploidy_priors} \
            --output ${output_dir_} \
            --output-prefix ${cohort_entity_id} \
            --verbosity DEBUG \
            --mean-bias-standard-deviation ${default="0.01" mean_bias_standard_deviation} \
            --mapping-error-rate ${default="0.01" mapping_error_rate} \
            --global-psi-scale ${default="0.001" global_psi_scale} \
            --sample-psi-scale ${default="0.0001" sample_psi_scale}

        tar czf ${cohort_entity_id}-contig-ploidy-model.tar.gz -C ${output_dir_}/${cohort_entity_id}-model .
        tar czf ${cohort_entity_id}-contig-ploidy-calls.tar.gz -C ${output_dir_}/${cohort_entity_id}-calls .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File contig_ploidy_model_tar = "${cohort_entity_id}-contig-ploidy-model.tar.gz"
        File contig_ploidy_calls_tar = "${cohort_entity_id}-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCohortMode {
    Int scatter_index
    String cohort_entity_id
    Array[File] read_count_files
    File contig_ploidy_calls_tar
    File intervals
    File? annotated_intervals
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    # Caller parameters
    Float? p_alt
    Float? p_active
    Float? cnv_coherence_length
    Float? class_coherence_length
    Int? max_copy_number

    # Denoising model parameters
    Int? max_bias_factors
    Float? mapping_error_rate
    Float? interval_psi_scale
    Float? sample_psi_scale
    Float? depth_correction_tau
    Float? log_mean_bias_standard_deviation
    Float? init_ard_rel_unexplained_variance
    Int? num_gc_bins
    Float? gc_curve_standard_deviation
    String? copy_number_posterior_expectation_mode
    Boolean? enable_bias_factors
    Int? active_class_padding_hybrid_mode

    # Hybrid ADVI parameters
    Float? learning_rate
    Float? adamax_beta_1
    Float? adamax_beta_2
    Int? log_emission_samples_per_round
    Float? log_emission_sampling_median_rel_error
    Int? log_emission_sampling_rounds
    Int? max_advi_iter_first_epoch
    Int? max_advi_iter_subsequent_epochs
    Int? min_training_epochs
    Int? max_training_epochs
    Float? initial_temperature
    Int? num_thermal_advi_iters
    Int? convergence_snr_averaging_window
    Float? convergence_snr_trigger_threshold
    Int? convergence_snr_countdown_window
    Int? max_calling_iters
    Float? caller_update_convergence_threshold
    Float? caller_internal_admixing_rate
    Float? caller_external_admixing_rate
    Boolean? disable_annealing

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=${default=8 cpu}
        export OMP_NUM_THREADS=${default=8 cpu}

        mkdir contig-ploidy-calls-dir
        tar xzf ${contig_ploidy_calls_tar} -C contig-ploidy-calls-dir

        gatk --java-options "-Xmx${command_mem_mb}m"  GermlineCNVCaller \
            --run-mode COHORT \
            -L ${intervals} \
            --input ${sep=" --input " read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls-dir \
            ${"--annotated-intervals " + annotated_intervals} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${output_dir_} \
            --output-prefix ${cohort_entity_id} \
            --verbosity DEBUG \
            --p-alt ${default="1e-6" p_alt} \
            --p-active ${default="1e-2" p_active} \
            --cnv-coherence-length ${default="10000.0" cnv_coherence_length} \
            --class-coherence-length ${default="10000.0" class_coherence_length} \
            --max-copy-number ${default="5" max_copy_number} \
            --max-bias-factors ${default="5" max_bias_factors} \
            --mapping-error-rate ${default="0.01" mapping_error_rate} \
            --interval-psi-scale ${default="0.001" interval_psi_scale} \
            --sample-psi-scale ${default="0.0001" sample_psi_scale} \
            --depth-correction-tau ${default="10000.0" depth_correction_tau} \
            --log-mean-bias-standard-deviation ${default="0.1" log_mean_bias_standard_deviation} \
            --init-ard-rel-unexplained-variance ${default="0.1" init_ard_rel_unexplained_variance} \
            --num-gc-bins ${default="20" num_gc_bins} \
            --gc-curve-standard-deviation ${default="1.0" gc_curve_standard_deviation} \
            --copy-number-posterior-expectation-mode ${default="HYBRID" copy_number_posterior_expectation_mode} \
            --enable-bias-factors ${default="true" enable_bias_factors} \
            --active-class-padding-hybrid-mode ${default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ${default="0.05" learning_rate} \
            --adamax-beta-1 ${default="0.9" adamax_beta_1} \
            --adamax-beta-2 ${default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ${default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ${default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ${default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ${default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ${default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ${default="10" min_training_epochs} \
            --max-training-epochs ${default="100" max_training_epochs} \
            --initial-temperature ${default="2.0" initial_temperature} \
            --num-thermal-advi-iters ${default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ${default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ${default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ${default="10" convergence_snr_countdown_window} \
            --max-calling-iters ${default="10" max_calling_iters} \
            --caller-update-convergence-threshold ${default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ${default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ${default="1.00" caller_external_admixing_rate} \
            --disable-annealing ${default="false" disable_annealing}

        tar czf ${cohort_entity_id}-gcnv-model-${scatter_index}.tar.gz -C ${output_dir_}/${cohort_entity_id}-model .
        CURRENT_SAMPLE=0
        while [ $CURRENT_SAMPLE -lt ${num_samples} ]; do
            tar czf ${cohort_entity_id}-shard-${scatter_index}-sample-$CURRENT_SAMPLE-gcnv-calls.tar.gz -C ${output_dir_}/${cohort_entity_id}-calls/SAMPLE_$CURRENT_SAMPLE .
            let CURRENT_SAMPLE=CURRENT_SAMPLE+1
        done
        tar czf ${cohort_entity_id}-gcnv-tracking-${scatter_index}.tar.gz -C ${output_dir_}/${cohort_entity_id}-tracking .
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File gcnv_model_tar = "${cohort_entity_id}-gcnv-model-${scatter_index}.tar.gz"
        File calling_config_json = "${output_dir_}/${cohort_entity_id}-calls/calling_config.json"
        File denoising_config_json = "${output_dir_}/${cohort_entity_id}-calls/denoising_config.json"
        File gcnvkernel_version_json = "${output_dir_}/${cohort_entity_id}-calls/gcnvkernel_version.json"
        File sharded_interval_list = "${output_dir_}/${cohort_entity_id}-calls/interval_list.tsv"
        Array[File] gcnv_call_tars = glob("*-gcnv-calls.tar.gz")
        File gcnv_tracking_tar = "${cohort_entity_id}-gcnv-tracking-${scatter_index}.tar.gz"
    }
}
