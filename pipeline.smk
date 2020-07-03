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
        fastqc --extract -o {output.folder} {input.reads[0]} {input.reads[1]} 1> {log.out_file} 2> {log.err_file}
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
        "bbmerge.sh in={input.filtered_reads} out={output.merged_reads} outu={output.unmerged_reads} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly__quick_assembly_with_tadpole"
rule assembly__quick_assembly_with_tadpole:
    # Static
    message:
        "Running step:" + rule_name
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
        contigs = rules.setup.params.folder + "/contigs.fasta"
    shell:
        "tadpole.sh in={input.merged_reads},{input.unmerged_reads} out={output.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__sketch_on_contigs"
rule assembly_check__sketch_on_contigs:
    # Static
    message:
        "Running step:" + rule_name
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__quick_assembly_with_tadpole.output
    output:
        sketch = rules.setup.params.folder + "/contigs.sketch"
    shell:
        "sendsketch.sh in={input.contigs} outsketch={output.sketch} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__stats"
rule post_assembly__stats:
    # Static
    message:
        "Running step:" + rule_name
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
        contigs = rules.assembly__quick_assembly_with_tadpole.output
    output:
        stats = touch(rules.setup.params.folder + "/post_assembly__stats")
    shell:
        "stats.sh {input.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__mapping_with_bbmap"
rule post_assembly__mapping_with_bbmap:
    # Static
    message:
        "Running step:" + rule_name
    shadow:
        "shallow"
    log:
        out_file = rules.setup.params.folder + "/log/" + rule_name + ".out.log",
        err_file = rules.setup.params.folder + "/log/" + rule_name + ".err.log",
    benchmark:
        rules.setup.params.folder + "/benchmarks/" + rule_name + ".benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__quick_assembly_with_tadpole.output,
        filtered = rules.setup__filter_reads_with_bbduk.output.filtered_reads
    output:
        mapped = temp(rules.setup.params.folder + "/contigs.sam")
    shell:
        "bbmap.sh ref={input.contigs} in={input.filtered} out={output.mapped} ambig=random 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__samtools_stats"
rule post_assembly__samtools_stats:
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
        mapped = rules.post_assembly__mapping_with_bbmap.output.mapped
    output:
        stats = rules.setup.params.folder + "/contigs.stats",
    shell:
        "samtools stats {input.mapped} 1> {output.stats} 2> {log.err_file}"


rule_name = "post_assembly__pileup"
rule post_assembly__pileup:
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
        mapped = rules.post_assembly__mapping_with_bbmap.output.mapped
    output:
        coverage = temp(rules.setup.params.folder + "/contigs.cov"),
        pileup = rules.setup.params.folder + "/contigs.pileup"
    shell:
        "pileup.sh in={input.mapped} basecov={output.coverage} out={output.pileup} 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__depth"
rule summarize__depth:
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
        coverage = rules.post_assembly__pileup.output.coverage
    params:
        sampleComponentObj = bifrost_sampleComponentObj
    output:
        contig_depth_yaml = rules.setup.params.folder + "/contigs.sum.cov",
        binned_depth_yaml = rules.setup.params.folder + "/contigs.bin.cov"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/rule__summarize_depth.py")


rule_name = "post_assembly__call_variants"
rule post_assembly__call_variants:
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
        contigs = rules.assembly__quick_assembly_with_tadpole.output,
        mapped = rules.post_assembly__mapping_with_bbmap.output.mapped
    output:
        variants = temp(rules.setup.params.folder + "/contigs.vcf"),
    shell:
        "callvariants.sh in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__variants"
rule summarize__variants:
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
        variants = rules.post_assembly__call_variants.output
    params:
        sampleComponentObj = bifrost_sampleComponentObj
    output:
        variants_yaml = rules.setup.params.folder + "/contigs.variants",
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "scripts/rule__summarize_variants.py")

rule_name = "rename_contigs"
rule rename_contigs:
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
        contigs = rules.assembly__quick_assembly_with_tadpole.output,
    output:
        contigs = rules.setup.params.folder + "/" + sample_name + ".fasta",
    params:
        sample_name = sample_name
    shell:
        "sed -e 's/Contig/{params.sample_name}/' {input.contigs} > {output.contigs}"
#* Dynamic section: end ****************************************************************************

rule_name = "datadump"
rule datadump:
    # Static
    message:
        "Running step:" + rule_name
    log:
        out_file = component_name + "/log/" + rule_name + ".out.log",
        err_file = component_name + "/log/" + rule_name + ".err.log",
    benchmark:
        component_name + "/benchmarks/" + rule_name + ".benchmark"
    input:
        #* Dynamic section: start ******************************************************************
        rules.fastqc_on_reads.output.fastqc_summary,
        rules.assembly_check__sketch_on_contigs.output.sketch,
        rules.post_assembly__stats.output.stats,
        rules.post_assembly__samtools_stats.output.stats,
        rules.summarize__depth.output.contig_depth_yaml,
        rules.summarize__variants.output.variants_yaml
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        sampleComponentObj = bifrost_sampleComponentObj
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
