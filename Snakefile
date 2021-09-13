configfile: "config.yaml"
bamPath = config["bamPath"]
outPath = config["outPath"]
cramPath = config["cramPath"]

SAMPLES = ["LP6008119-DNA_B04","LP6008119-DNA_E11"]

rule CramToBam:
    input:
        cram_file=cramPath + "{sample}.cram"
    output:
        bam_file=bamPath + "{sample}.bam",
        bam_index=bamPath + "{sample}.bam.bai",
    benchmark:
        "benchmarks/CramToBam/{sample}.tsv"
    conda:
        "envs/sam_only.yaml"
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        samtools view -b -h -@ {threads} -T {config[hg19ref]} -o {output.bam_file} {input.cram_file}
        samtools index -@ {threads} {output.bam_file}
        """

rule picard:
    input:
        bamPath + "{sample}.bam"
    output:
        outPath + "sorted/{sample}.bam"
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
        outPath + "sam/{sample}.sam"
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
    envmodules:
        "apps/steak/20190912-singularity"
    shell:
        "steak --input {input} --TE-reference  {config[HERVK_ref]} --paired --aligned --output {config[outPath]}steak/{wildcards.sample}.Steak.txt"

rule novoalign_trimmed:
    input:
       outPath + "steak/{sample}.Steak.txt.te.fastq"
    output:
       outPath + "novoalign_trimmed/{sample}.bam"
    threads: 8
    conda:
      "envs/novo.yaml"
    benchmark:
        "benchmarks/{sample}.novo_trimmed.benchmark.txt"
    log:
       "logs/novoalign_trimmed/{sample}.log"
    shell:
       """
       if [ -f {config[novoindex]} ]; then
          echo "novoalign index found"
       else 
           novoindex {config[novoindex]} {config[hg19ref]} 
       fi
       novoalign -d {config[novoindex]} -o SAM -f {input} | samtools view -Sb  > {output} 
       """ 

rule novoalign_guided:
    input:
       outPath + "steak/{sample}.Steak.txt.1.fastq",
       outPath + "steak/{sample}.Steak.txt.2.fastq"
      # outPath + "steak/{sample}.Steak.txt" #edit that, it has 1.fastq etc appended
    output:
        outPath + "novoalign_guided/{sample}.bam"
    threads: 8
    benchmark:
        "benchmarks/{sample}.novo_guided.benchmark.txt"
    conda:
      "envs/novo.yaml"
    log:
       "logs/novoalign_guided/{sample}.log"
    shell:
       """
       if [ -f {config[novoindex]} ]; then
           echo "novoalign index found"
       else 
           novoindex {config[novoindex]} {config[hg19ref]} 
       fi
       novoalign -d {config[novoindex]} -o SAM -o FullNW -f {input} | samtools view -Sb > {output}
       """




rule bamToBedTrimmed:
     input:
      outPath + "novoalign_trimmed/{sample}.bam"
     output:
      outPath + "bamToBed_trimmed/{sample}.bed"
     log:
       "logs/bamToBed_trimmed/{sample}.log"
     conda:
      "envs/bedtools.yaml"
     shell:
       "bedtools bamtobed -i {input} > {output}"

rule bamToBedGuided:
     input:
      outPath + "novoalign_guided/{sample}.bam"
     output:
      outPath + "bamToBed_guided/{sample}.bed"
     log:
       "logs/bamToBed_guided/{sample}.log"
     conda:
      "envs/bedtools.yaml"
     shell:
       "bedtools bamtobed -i {input} > {output}"


rule markKnown:
     input:
        G=outPath + "bamToBed_guided/{sample}.bed",
        T=outPath + "bamToBed_trimmed/{sample}.bed"
     output:
        outPath + "known/{sample}.knownHits.bed"
     log:
        "logs/known/{sample}.log"
     conda:
      "envs/bedtools.yaml"
     shell:
        """
        bedtools window -w 100 -c -a {config[knownNR]} -b  {input.G}  > {wildcards.sample}.guided.readCounts.bed 
        bedtools window -w 100 -c -a {config[knownNR]} -b  {input.T}  > {wildcards.sample}.trimmed.readCounts.bed 
        python {config[pythonScripts]}/mergeTrimmedAndGuided.py {wildcards.sample}.guided.readCounts.bed {wildcards.sample}.trimmed.readCounts.bed > {output} 
        """

rule markNovel:
     input:
        G=outPath + "bamToBed_guided/{sample}.bed",
      	T=outPath + "bamToBed_trimmed/{sample}.bed"
     output:
        outPath + "novel/{sample}.novelHits.bed"
     log:
        "logs/novel/{sample}.log"
     conda:
      "envs/bedtools.yaml"
     shell:
        """
        bedtools sort -i  {input.T} >  {config[outPath]}novel/{wildcards.sample}.sortedTrimmed.bed 
        bedtools window -w 2000 -v -a {config[outPath]}novel/{wildcards.sample}.sortedTrimmed.bed -b {config[knownNR]}| bedtools merge -c 4 -o count -d 1000 -i - | awk '{{if ($4 >=5) print ($0);}}'  > {config[outPath]}novel/{wildcards.sample}.trimmed_novelFiltered.bed 
        bedtools sort -i  {input.G} >  {config[outPath]}novel/{wildcards.sample}.sortedGuided.bed 
        bedtools window -w 2000 -v -a {config[outPath]}novel/{wildcards.sample}.sortedGuided.bed  -b {config[knownNR]}| bedtools merge -c 4 -o count -d 1000 -i - | awk '{{if ($4 >=5) print ($0);}}'  > {config[outPath]}novel/{wildcards.sample}.guided_novelFiltered.bed 
        cat {config[outPath]}novel/{wildcards.sample}.guided_novelFiltered.bed {config[outPath]}novel/{wildcards.sample}.trimmed_novelFiltered.bed  > {config[outPath]}novel/{wildcards.sample}.novelFiltered.bed 
        bedtools sort -i {config[outPath]}novel/{wildcards.sample}.novelFiltered.bed   >  {config[outPath]}novel/{wildcards.sample}.novelFiltered.sorted.bed 
        rm {config[outPath]}novel/{wildcards.sample}.novelFiltered.bed 
        bedtools merge -i {config[outPath]}novel/{wildcards.sample}.novelFiltered.sorted.bed > {output} 
        rm {config[outPath]}novel/{wildcards.sample}.novelFiltered.sorted.bed 
        """

rule merge:
    input:
        bed=expand(outPath + "bamToBed_guided/{sample}.bed", sample=SAMPLES)
    output:
        outPath + "final/mergedBed.bed"
    conda:
      "envs/bedtools.yaml"
    shell:
        #"cat {input.bed} > {output}"
        "cat {input.bed} | grep -E 'chr[[:digit:]XY]' > {config[outPath]}final/noHeader.bed "
        "bedtools sort -i {config[outPath]}final/noHeader.bed > {config[outPath]}final/noHeader.sorted.bed"
        "bedtools merge -i {config[outPath]}final/noHeader.sorted.bed > {output}"
        "rm {config[outPath]}final/noHeader.bed"
        "rm {config[outPath]}final/noHeader.sorted.bed"
