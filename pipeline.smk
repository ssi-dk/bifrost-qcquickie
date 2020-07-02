import os
from bifrostlib import datahandling


os.umask(0o2)
bifrost_sampleComponentObj = datahandling.SampleComponentObj(config["sample_id"], config["component_id"])
sample_name, component_name, dockerfile, options, bifrost_resources = bifrost_sampleComponentObj.load()
bifrost_sampleComponentObj.started()


onerror:
   bifrost_sampleComponentObj.failure()


rule all:
    input:
         # file is defined by datadump function
        component_name + "/datadump_complete"


rule setup:
    output:
        init_file = touch(
            temp(component_name + "/initialized")),
    params:
        folder = component_name


rule_name = "check_requirements"
rule check_requirements:
    message:
        "Running step:" + rule_name
    log:
        out_file = component_name + "/log/" + rule_name + ".out.log",
        err_file = component_name + "/log/" + rule_name + ".err.log",
    benchmark:
        component_name + "/benchmarks/" + rule_name + ".benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = component_name + "/requirements_met",
    params:
        bifrost_sampleComponentObj
    run:
        bifrost_sampleComponentObj.check_requirements()
#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************

rule_name = "fastqc_on_reads"
rule fastqc_on_reads:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.check_requirements.output.check_file,
        reads = bifrost_sampleComponentObj.get_reads()
    output:
        folder = directory(rules.setup.params.folder + "/fastqc"),
        fastqc_summary = rules.setup.params.folder + "/fastqc_data.txt"
    shell:
        """
        mkdir {output.folder}
        fastqc --extract -o {output.folder} -t {threads} {input.reads[0]} {input.reads[1]} 1> {log.out_file} 2> {log.err_file}
        cat {output.folder}/*/fastqc_data.txt > {output.fastqc_summary}
        """


rule_name = "setup__filter_reads_with_bbduk"
rule setup__filter_reads_with_bbduk:
    # Static
    message:
        "Running step:" + rule_name
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        folder = rules.check_requirements.output.check_file,
        reads = bifrost_sampleComponentObj.get_reads()
    output:
        filtered_reads = temp(rules.setup.params.folder + "/filtered.fastq")
    params:
        adapters = bifrost_resources["adapters_fasta"]
    shell:
        "bbduk.sh in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo qtrim=r minlength=30 json=t 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__combine_reads_with_bbmerge"
rule assembly_check__combine_reads_with_bbmerge:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        filtered_reads = rules.setup__filter_reads_with_bbduk.output.filtered_reads,
    output:
        merged_reads = temp(rules.setup.params.folder + "/merged.fastq"),
        unmerged_reads = temp(rules.setup.params.folder + "/unmerged.fastq")
    shell:
        "bbmerge.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.filtered_reads} out={output.merged_reads} outu={output.unmerged_reads} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__quick_assembly_with_tadpole"
rule assembly_check__quick_assembly_with_tadpole:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        merged_reads = rules.assembly_check__combine_reads_with_bbmerge.output.merged_reads,
        unmerged_reads = rules.assembly_check__combine_reads_with_bbmerge.output.unmerged_reads
    output:
        contigs = temp(rules.setup.params.folder + "/raw_contigs.fasta")
    shell:
        "tadpole.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.merged_reads},{input.unmerged_reads} out={output.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__rename_contigs"
rule assembly_check__rename_contigs:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.assembly_check__quick_assembly_with_tadpole.output.contigs,
    output:
        contigs = rules.setup.params.folder + "/contigs.fasta"
    run:
        try:
            input_file = str(input.contigs)
            output_file = str(output.contigs)
            log_out = str(log.out_file)
            log_err = str(log.err_file)

            datahandling.log(log_out, "Started {}\n".format(rule_name))
            with open(input_file, "r") as fasta_input:
                records = list(Bio.SeqIO.parse(fasta_input, "fasta"))
            for record in records:
                record.id = record.id.split(",")[0]
                record.description = record.id
            with open(output_file, "w") as output_handle:
                Bio.SeqIO.write(records, output_handle, "fasta")
            datahandling.log(log_out, "Done {}\n".format(rule_name))
        except Exception as e:
            datahandling.log(log_err, str(traceback.format_exc()))


rule_name = "assembly_check__sketch_on_contigs"
rule assembly_check__sketch_on_contigs:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.assembly_check__rename_contigs.output.contigs
    output:
        sketch = rules.setup.params.folder + "/contigs.sketch"
    shell:
        "sendsketch.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.contigs} outsketch={output.sketch} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__stats"
rule post_assembly__stats:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    message:
        "Running step: {rule}"
    input:
        contigs = rules.assembly_check__rename_contigs.output.contigs
    output:
        stats = touch(rules.setup.params.folder + "/post_assermbly__stats")
    shell:
        "stats.sh -Xmx{resources.memory_in_GB}G {input.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__map_reads_to_assembly_with_bbmap"
rule assembly_check__map_reads_to_assembly_with_bbmap:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.assembly_check__rename_contigs.output.contigs,
        filtered = rules.setup__filter_reads_with_bbduk.output.filtered_reads
    output:
        mapped = temp(rules.setup.params.folder + "/contigs.sam")
    shell:
        "bbmap.sh threads={threads} -Xmx{resources.memory_in_GB}G ref={input.contigs} in={input.filtered} out={output.mapped} ambig=random 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__samtools_stats"
rule post_assembly__samtools_stats:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        mapped = rules.assembly_check__map_reads_to_assembly_with_bbmap.output.mapped
    output:
        stats = rules.setup.params.folder + "/contigs.stats",
    shell:
        "samtools stats -@ {threads} {input.mapped} 1> {output.stats} 2> {log.err_file}"


rule_name = "assembly_check__pileup_on_mapped_reads"
rule assembly_check__pileup_on_mapped_reads:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        mapped = rules.assembly_check__map_reads_to_assembly_with_bbmap.output.mapped
    output:
        coverage = temp(rules.setup.params.folder + "/contigs.cov"),
        pileup = rules.setup.params.folder + "/contigs.pileup"
    shell:
        "pileup.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.mapped} basecov={output.coverage} out={output.pileup} 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__depth"
rule summarize__depth:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        coverage = rules.assembly_check__pileup_on_mapped_reads.output.coverage
    output:
        contig_depth_yaml = rules.setup.params.folder + "/contigs.sum.cov",
        binned_depth_yaml = rules.setup.params.folder + "/contigs.bin.cov"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/summarize_depth.py")


rule_name = "assembly_check__call_variants"
rule assembly_check__call_variants:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.assembly_check__rename_contigs.output.contigs,
        mapped = rules.assembly_check__map_reads_to_assembly_with_bbmap.output.mapped,
    output:
        variants = temp(rules.setup.params.folder + "/contigs.vcf")
    shell:
        "callvariants.sh threads={threads} -Xmx{resources.memory_in_GB}G in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__variants"
rule summarize__variants:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        variants = rules.assembly_check__call_variants.output.variants,
    output:
        variants_yaml = rules.setup.params.folder + "/contigs.variants",
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/summarize_variants.py")


rule_name = "datadump_qcquickie"
rule datadump_qcquickie:
    # Static
    message:
        "Running step:" + rule_name
    threads:
        global_threads
    resources:
        memory_in_GB = global_memory_in_GB
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        rules.fastqc_on_reads.output.fastqc_summary,
        rules.assembly_check__sketch_on_contigs.output.sketch,
        rules.post_assembly__stats.output.stats,
        rules.post_assembly__samtools_stats.output.stats,
        rules.summarize__depth.output.contig_depth_yaml,
        rules.summarize__variants.output.variants_yaml,
    output:
        summary = touch(rules.all.input)
    params:
        sample = db_sample.get("name", "ERROR") + "__" + component + ".yaml",
        folder = rules.setup.params.folder,
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
