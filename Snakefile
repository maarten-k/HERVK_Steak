configfile: "config.yaml"


containerized: config["condacontainer"]


outPath = config["outPath"]
cramPath = config["cramPath"]

SAMPLES = [config["sample"]]


rule all:
    input:
        expand(outPath + "bamToBed/{sample}.guided.bed", sample=SAMPLES),
        expand(outPath + "bamToBed/{sample}.trimmed.bed", sample=SAMPLES),
        expand(outPath + "known/{sample}.knownHits.bed", sample=SAMPLES),
        expand(outPath + "novel/{sample}.novelHits.bed", sample=SAMPLES),


rule SortOnName:
    input:
        cramPath + "/{sample}.cram",
    output:
        temp(outPath + "sam/{sample}.sqfs"),
    group:"squash"
    threads: 1
    benchmark:
        #repeat("benchmarks/{sample}.picard.benchmark.txt",3)  
        "benchmarks/{sample}.picard.benchmark.txt"
    conda:
        "envs/sortonname.yaml"
    log:
        "logs/picard/{sample}.log",
    shell:
        """
        picard SortSam -I {input} -O /dev/stdout -REFERENCE_SEQUENCE hg19.fa -COMPRESSION_LEVEL 0 -SORT_ORDER queryname -VALIDATION_STRINGENCY LENIENT -QUIET true | samtools view -O sam --input-fmt-option required_fields=0x23d |mksquashfs /dev/null {output} -comp zstd -Xcompression-level 1 -p "{wildcards.sample}.sam f 644 0 0 cat" 
        """


rule MountSquasfs:
    input:
        outPath + "sam/{sample}.sqfs",
    output:
        temp(outPath + "sam/{sample}_mounted"),
    params:
        sam=outPath + "sam/{sample}_mount/{sample}.sam",
        dir=directory(outPath + "sam/{sample}_mount"),
    group:"squash"
    conda:
        "envs/squashfs.yaml"
    shell:
        """
        mkdir -p {params.dir}
        squashfuse {input} {params.dir}
        touch {output} 
        """


rule steak:
    input:
        sqfs=outPath + "sam/{sample}.sqfs",
        mounted=outPath + "sam/{sample}_mounted",
        steaksif=config["steakcontainer"],
    output:
        outPath + "steak/{sample}.Steak.txt.1.fastq",
        outPath + "steak/{sample}.Steak.txt.2.fastq",
        outPath + "steak/{sample}.Steak.txt.te.fastq",
    params:
        sam=outPath + "sam/{sample}_mount/{sample}.sam",
    group:"squash"
    threads: 1
    benchmark:
        "benchmarks/{sample}.steak.benchmark.txt"
    log:
        "logs/steak/{sample}.log",
    container:
        None
    shell:
        """
        /usr/bin/singularity exec {input.steaksif}  /STEAK-master/steak  --input {params.sam}  --TE-reference  {config[HERVK_ref]} --paired --aligned --output {config[outPath]}steak/{wildcards.sample}.Steak.txt
        """


rule unmountSqaushfs:
    input:
        outPath + "steak/{sample}.Steak.txt.1.fastq",
        lockfile=outPath + "sam/{sample}_mounted",
        sqfs=outPath + "sam/{sample}.sqfs",
    params:
        dir=outPath + "sam/{sample}_mount",
    group:"squash"
    output:
        temp(outPath + "sam/{sample}.unmount"),
    conda:
        "envs/squashfs.yaml"
    shell:
        """
        fusermount -u {params.dir}
        #rm {input.lockfile}
        touch {output}
        """


rule novoalign:
    input:
        A=outPath + "steak/{sample}.Steak.txt.1.fastq",
        B=outPath + "steak/{sample}.Steak.txt.2.fastq",
        C=outPath + "steak/{sample}.Steak.txt.te.fastq",
        unmount=outPath + "sam/{sample}.unmount",
    output:
        outPath + "novoalign/{sample}.guided.bam",
        outPath + "novoalign/{sample}.trimmed.bam",
    threads: 1
    conda:
        "envs/novo.yaml"
    benchmark:
        "benchmarks/{sample}.novoalign.benchmark.txt"
    log:
        "logs/novoalign/{sample}.log",
    shell:
        """
        novoalign -d {config[novoindex]} -o SAM -o FullNW -f {input.A} {input.B} | samtools view -Sb > {output[0]}
        novoalign -d {config[novoindex]} -o SAM -f {input.C} | samtools view -Sb  > {output[1]}
        """


rule bamToBed:
    input:
        outPath + "novoalign/{sample}.guided.bam",
        outPath + "novoalign/{sample}.trimmed.bam",
    output:
        outPath + "bamToBed/{sample}.guided.bed",
        outPath + "bamToBed/{sample}.trimmed.bed",
    log:
        "logs/bamToBed/{sample}.log",
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools bamtobed -i {input[0]} > {output[0]} 
        bedtools bamtobed -i {input[1]} > {output[1]} 
        """


rule markKnown:
    input:
        G=outPath + "bamToBed/{sample}.guided.bed",
        T=outPath + "bamToBed/{sample}.trimmed.bed",
    output:
        outPath + "known/{sample}.knownHits.bed",
    log:
        "logs/known/{sample}.log",
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools window -w 500 -c -a {config[knownNR]} -b  {input.G}  > {config[outPath]}known/{wildcards.sample}.guided.readCounts.bed 
        bedtools window -w 500 -c -a {config[knownNR]} -b  {input.T}  >  {config[outPath]}known/{wildcards.sample}.trimmed.readCounts.bed 
        python {config[pythonScripts]}/mergeTrimmedAndGuided.py  {config[outPath]}known/{wildcards.sample}.guided.readCounts.bed  {config[outPath]}known/{wildcards.sample}.trimmed.readCounts.bed > {output} 
        rm  {config[outPath]}known/{wildcards.sample}.guided.readCounts.bed  
        rm  {config[outPath]}known/{wildcards.sample}.trimmed.readCounts.bed  
        """


rule markNovel:
    input:
        G=outPath + "bamToBed/{sample}.guided.bed",
        T=outPath + "bamToBed/{sample}.trimmed.bed",
    output:
        outPath + "novel/{sample}.novelHits.bed",
    log:
        "logs/novel/{sample}.log",
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools sort -i  {input.T} >  {config[outPath]}novel/{wildcards.sample}.sortedTrimmed.bed 
        bedtools window -w 500 -v -a {config[outPath]}novel/{wildcards.sample}.sortedTrimmed.bed -b {config[knownNR]}| bedtools merge -c 4 -o count -d 1000 -i - | awk '{{if ($4 >=5) print ($0);}}'  > {config[outPath]}novel/{wildcards.sample}.trimmed_novelFiltered.bed 
        bedtools sort -i  {input.G} >  {config[outPath]}novel/{wildcards.sample}.sortedGuided.bed 
        bedtools window -w 500 -v -a {config[outPath]}novel/{wildcards.sample}.sortedGuided.bed  -b {config[knownNR]}| bedtools merge -c 4 -o count -d 1000 -i - | awk '{{if ($4 >=5) print ($0);}}'  > {config[outPath]}novel/{wildcards.sample}.guided_novelFiltered.bed 
        cat {config[outPath]}novel/{wildcards.sample}.guided_novelFiltered.bed {config[outPath]}novel/{wildcards.sample}.trimmed_novelFiltered.bed  > {config[outPath]}novel/{wildcards.sample}.novelFiltered.bed 
        bedtools sort -i {config[outPath]}novel/{wildcards.sample}.novelFiltered.bed   >  {config[outPath]}novel/{wildcards.sample}.novelFiltered.sorted.bed 
        rm {config[outPath]}novel/{wildcards.sample}.novelFiltered.bed 
        bedtools merge -i {config[outPath]}novel/{wildcards.sample}.novelFiltered.sorted.bed > {output} 
        rm {config[outPath]}novel/{wildcards.sample}.novelFiltered.sorted.bed 
        rm {config[outPath]}novel/{wildcards.sample}.sortedTrimmed.bed 
        rm {config[outPath]}novel/{wildcards.sample}.sortedGuided.bed 
        rm {config[outPath]}novel/{wildcards.sample}.trimmed_novelFiltered.bed 
        rm {config[outPath]}novel/{wildcards.sample}.guided_novelFiltered.bed  
        """
