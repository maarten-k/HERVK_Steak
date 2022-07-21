configfile: "config.yaml"
bamPath = config["bamPath"]
outPath = config["outPath"]
cramPath = config["cramPath"]

SAMPLES = ["samplecram"]

rule all:
    input:
        expand(outPath + "bamToBed/{sample}.guided.bed", sample=SAMPLES),
        expand(outPath + "bamToBed/{sample}.trimmed.bed", sample=SAMPLES),
        expand(outPath + "known/{sample}.knownHits.bed",sample=SAMPLES),
	expand(outPath + "novel/{sample}.novelHits.bed",sample=SAMPLES)

rule CramToBam:
    input:
        cram_file=cramPath + "{sample}.cram"
    output:
        temp(bamPath + "{sample}.bam"), 
        bamPath + "{sample}.bam.bai"
    benchmark:
        "benchmarks/CramToBam/{sample}.tsv"
    conda:
        "envs/sam_only.yaml"
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        samtools view -b -h -@ {threads} -T {config[hg19ref]} -o {output[0]} {input.cram_file}
        samtools index -@ {threads} {output[0]}
        """

rule picard:
    input:
        cramPath + "{sample}.cram"
    output:
        temp(outPath + "sorted/{sample}.bam")
    threads: 8
    benchmark:
      #repeat("benchmarks/{sample}.picard.benchmark.txt",3)  
      "benchmarks/{sample}.picard.benchmark.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/{sample}.log"
    shell:
        "picard SortSam I={input} O={output} SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT > {log}"

rule bamToSam:
    input:
        outPath + "sorted/{sample}.bam"
    output:
        temp(outPath + "sam/{sample}.sam")
    threads: 8
    conda:
        "envs/sam_only.yaml"
    log:
       "logs/sam/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bamToSam.benchmark.txt"
        #repeat("benchmarks/{sample}.bamToSam.benchmark.txt",3)
    shell:
        "samtools view {input} -h -o {output}"

rule steak:
    input:
       outPath + "sam/{sample}.sam"
    output:
       outPath + "steak/{sample}.Steak.txt.1.fastq",
       outPath + "steak/{sample}.Steak.txt.2.fastq",
       outPath + "steak/{sample}.Steak.txt.te.fastq"
    threads: 8
    benchmark:
        "benchmarks/{sample}.steak.benchmark.txt"
        #repeat("benchmarks/{sample}.steak.benchmark.txt",3)
    log:
       "logs/steak/{sample}.log"
    shell:
        "singularity exec ../steak.sif /STEAK-master/steak  --input {input} --TE-reference  {config[HERVK_ref]} --paired --aligned --output {config[outPath]}steak/{wildcards.sample}.Steak.txt"

rule novoalign:
    input:
       A=outPath + "steak/{sample}.Steak.txt.1.fastq",
       B=outPath + "steak/{sample}.Steak.txt.2.fastq",
       C=outPath + "steak/{sample}.Steak.txt.te.fastq"
    output:
       outPath + "novoalign/{sample}.guided.bam",
       outPath + "novoalign/{sample}.trimmed.bam"
    threads: 8
    conda:
      "envs/novo.yaml"
    benchmark:
        "benchmarks/{sample}.novoalign.benchmark.txt"
    log:
       "logs/novoalign/{sample}.log"
    shell:
       """
       novoalign -d {config[novoindex]} -o SAM -o FullNW -f {input.A} {input.B} | samtools view -Sb > {output[0]}
       novoalign -d {config[novoindex]} -o SAM -f {input.C} | samtools view -Sb  > {output[1]}
       """

rule bamToBed:
     input:
      outPath + "novoalign/{sample}.guided.bam",
      outPath + "novoalign/{sample}.trimmed.bam"
     output:
      outPath + "bamToBed/{sample}.guided.bed",
      outPath + "bamToBed/{sample}.trimmed.bed"
     log:
       "logs/bamToBed/{sample}.log"
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
        T=outPath + "bamToBed/{sample}.trimmed.bed"
     output:
        outPath + "known/{sample}.knownHits.bed"
     log:
        "logs/known/{sample}.log"
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
      	T=outPath + "bamToBed/{sample}.trimmed.bed"
     output:
        outPath + "novel/{sample}.novelHits.bed"
     log:
        "logs/novel/{sample}.log"
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
