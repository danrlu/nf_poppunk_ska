#!/usr/bin/env nextflow 

nextflow.enable.dsl=2

params.input_folder = "${workflow.launchDir}"
params.sample_sheet = "${params.input_folder}/r_file.txt"
params.plot_rmd = "${workflow.projectDir}/plot.Rmd"
params.ouput_folder = "${params.input_folder}/pop_ska"

plot_rmd = Channel.fromPath(params.plot_rmd, checkIfExists: true)

sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)

fq_in = sample_sheet.splitText()

/*
if (params.se) {

  fq_in = sample_sheet.splitCsv(header:false, sep: "\t")
                      .map { row -> [row[0], file(row[1]), file(row[2])]}

} else {

  fq_in = sample_sheet.splitCsv(header:false, sep: "\t")
                      .map { row -> [row[0], file(row[1])]}

}
*/


workflow { 

    sample_sheet | poppunk

    fq_in | ska_sketch 

    ska_sketch.out.collect() | ska_summary

//    poppunk.out.combine(ska_summary.out).map{[it]}.combine(plot_rmd) | plot

    poppunk.out.combine(ska_summary.out).combine(plot_rmd) | plot
}



process poppunk {
//    debug true

//    memory 16.GB

    publishDir "${params.ouput_folder}", mode: 'copy', pattern: "*"

    input:
      file(ss)

    output:
      path "dist.txt", emit: dist
//      path "*.png"

    """
#    head $ss
    sketchlib sketch -l $ss -o database -s 10000 -k 15,31,2 --cpus 4
    sketchlib query dist database --cpus 4 > dist.txt

    """
}



process ska_sketch {

//    debug true


//    publishDir "${workflow.projectDir}/${params.input_folder}/output_ska", mode: 'copy', pattern: "*"

    input:
      val(fq_in)

    output:
      path "*.skf"
//      path "*.png"

    shell:
    '''
    s=($(echo -n !{fq_in}))
#    echo $s
#    echo ${#s[@]}
#    echo ${s[0]} ${s[1]}
    if [ "${#s[@]}" -eq 2 ]; then
      ska fastq -o ${s[0]} ${s[1]}
    elif [ "${#s[@]}" -eq 3 ]; then
      ska fastq -o ${s[0]} ${s[1]} ${s[2]}
    fi

    '''
}


process ska_summary {

    publishDir "${params.ouput_folder}", mode: 'copy', pattern: "*"

    input:
      path("*")

    output:
      path "*"

    """
    ska summary `ls *.skf` > summary.txt

    ska distance -s 25 -i 0.95 `ls *.skf` > distance.tsv

    dot -Tpng distances.dot -o distances.png
    """

}



process plot {

//    publishDir "${params.input_folder}/output_pp", mode: 'copy', pattern: '*.html'
    conda "/home/genepi-powerusers/miniconda3/envs/r"

    input:
        path("*")
//      tuple path("*"), path(rmd)

    """
    cat plot.Rmd > ${params.ouput_folder}/plot.Rmd
#    cp -H plot.Rmd ${params.ouput_folder}/plot.Rmd
    cd ${params.ouput_folder}
#    Rscript -e "rmarkdown::render('plot.Rmd')"
    Rscript -e "rmarkdown::render('plot.Rmd')"
    """
}
