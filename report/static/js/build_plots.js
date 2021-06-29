// check if elemennt exists in array
function checkAvailability(arr, val) {
    return arr.some(function (arrVal) {
        return val === arrVal;
    });
}
// build all plots
function build_plot(){
    //data quality
    build_mapping_plot();
    build_coverage_plot();
    build_mean_quality_plot();
    build_gc_content_plot();
    build_insert_size_plot();
    //copy number variation
    build_clonality_plot();
    build_purity_plot();
    build_ploidy_plot();
    //somatic variants
    build_sv_summary_plots();
}

function build_sv_summary_plots(){
    build_variant_classification_plot();
    build_variant_type_plot();
    build_snv_class_plot();
    build_variant_per_sample_plot();
    build_variant_classification_summary_plot();
    build_top_10_genes_plot();
}
// variant classification types
var vc_types = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Splice_Site','In_Frame_Del','Frame_Shift_Ins','In_Frame_Ins','Translation_Start_Site','Nonstop_Mutation'];
var bases = ['T','C','G','A']
// snv class conversion from maftools
var snv_conv = {
    'A>G':'T>C',
    'T>C':'T>C',
    'C>T':'C>T',
    'G>A':'C>T',
    'A>T':'T>A',
    'T>A':'T>A',
    'A>C':'T>G',
    'T>G':'T>G',
    'C>A':'C>A',
    'G>T':'C>A',
    'C>G':'C>G',
    'G>C':'C>G'
  };
// find indices of top values
function findIndicesOfMax(inp, count) {
    var outp = new Array();
    for (var i = 0; i < inp.length; i++) {
        outp.push(i);
        if (outp.length > count) {
            outp.sort(function(a, b) { return inp[b] - inp[a]; });
            outp.pop();
        }
    }
    return outp;
}

// Data quality plots

var sample_type = ['normal', 'tumor'];
var sample_type_suffix = ['.N', '.T'];

function build_mapping_plot() {

    let total_reads = [];
    let mapped_reads = [];
    let dedup_reads = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // different types of reads - counts
                    total_reads.push(wes_data[i][sample_type[j]]['alignment']['total_reads'])
                    mapped_reads.push(wes_data[i][sample_type[j]]['alignment']['mapped_reads'])
                    dedup_reads.push(wes_data[i][sample_type[j]]['alignment']['dedup_reads'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }

    let trace1 = {
        x: sample_ids,
        y: total_reads,
        name: 'Total Reads',
        type: 'bar'
    };

    let trace2 = {
        x: sample_ids,
        y: mapped_reads,
        name: 'Mapped Reads',
        type: 'bar'
    };

    let trace3 = {
        x: sample_ids,
        y: dedup_reads,
        name: 'Dedup Reads',
        type: 'bar'
    };

    let data = [trace1, trace2, trace3];

    let layout = {
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Reads' },
            automargin: true
        },
    };

    Plotly.newPlot("mapping_plot", data, layout);

};

function build_coverage_plot() {

    let mean_depth = [];
    let sample_ids = [];

    let hover_text = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // read depth coverage
                    mean_depth.push(wes_data[i][sample_type[j]]['coverage']['mean_depth'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                    // extra information for hover
                    hover_text.push(
                        'Sample ID: ' + wes_data[i]['id'] + sample_type_suffix[j] + '<br>' +
                        'Mean Depth: ' + wes_data[i][sample_type[j]]['coverage']['mean_depth'] + '<br>' +
                        'Median Depth: ' + wes_data[i][sample_type[j]]['coverage']['median_depth'] + '<br>' +
                        'Q1 Depth: ' + wes_data[i][sample_type[j]]['coverage']['q1_depth'] + '<br>' +
                        'Q3 Depth: ' + wes_data[i][sample_type[j]]['coverage']['q3_depth'] + '<br>' +
                        '% Bases > 50: ' + wes_data[i][sample_type[j]]['coverage']['percent_bases_gt_50']
                    )
                }
            }
        }
    }

    let trace1 = {
        y: sample_ids,
        x: mean_depth,
        name: 'Mean Depth',
        type: 'bar',
        orientation: 'h',
        hovertemplate: '%{text}<extra></extra>',
        text: hover_text,
    };

    let data = [trace1];

    let layout = {
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Mean Depth' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
    };

    Plotly.newPlot("coverage_plot", data, layout);

};

function build_gc_content_plot() {

    let gc_content = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // gc content
                    gc_content.push({
                        y: wes_data[i][sample_type[j]]['alignment']['gc_content'],
                        mode: 'lines',
                        type: 'scatter',
                        name: wes_data[i]['id'] + sample_type_suffix[j]
                    })
                }
            }
        }
    }

    let data = gc_content;

    let layout = {
        barmode: 'overlay',
        hovermode: 'closest',
        xaxis: {
            title: { text: '% GC bases' },
            automargin: true
        },
        yaxis: {
            title: { text: 'GC Content' },
            automargin: true
        },
        legend: { title: { "text": "Sample" } }
    };

    Plotly.newPlot("gc_content_plot", data, layout);

};

function build_insert_size_plot() {

    let insert_size = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // insert size
                    insert_size.push({
                        y: wes_data[i][sample_type[j]]['alignment']['insert_size'],
                        mode: 'lines',
                        type: 'scatter',
                        name: wes_data[i]['id'] + sample_type_suffix[j]
                    })
                }
            }
        }
    }

    let data = insert_size;

    let layout = {
        barmode: 'overlay',
        hovermode: 'closest',
        xaxis: {
            title: { text: 'Insert Size' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Counts' },
            automargin: true
        },
        legend: { title: { "text": "Sample" } }
    };

    Plotly.newPlot("insert_size_plot", data, layout);

};



function build_mean_quality_plot() {

    let mean_quality = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // mean quality score
                    mean_quality.push(wes_data[i][sample_type[j]]['alignment']['mean_quality_score'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }

    let trace1 = {
        x: sample_ids,
        y: mean_quality,
        type: 'bar'
    };

    let data = [trace1];

    let layout = {
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Mean Quality Score' },
            automargin: true
        },
    };

    Plotly.newPlot("mean_quality_plot", data, layout);

};

// CNV plots

function build_clonality_plot() {

    let clonality = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                clonality.push(wes_data[i]['copy_number']['clonality'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    let trace1 = {
        x: sample_ids,
        y: clonality,
        type: 'bar'
    };

    let data = [trace1];

    let layout = {
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Clonality' },
            automargin: true
        },
    };

    Plotly.newPlot("clonality_plot", data, layout);

};

function build_purity_plot() {

    let purity = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                purity.push(wes_data[i]['copy_number']['purity'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    let trace1 = {
        x: sample_ids,
        y: purity,
        type: 'bar'
    };

    let data = [trace1];

    let layout = {
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Purity' },
            automargin: true
        },
    };

    Plotly.newPlot("purity_plot", data, layout);

};

function build_ploidy_plot() {

    let ploidy = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                ploidy.push(wes_data[i]['copy_number']['ploidy'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    let trace1 = {
        x: sample_ids,
        y: ploidy,
        mode: 'lines',
        type: 'scatter',
    };

    let data = [trace1];

    let layout = {
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Ploidy' },
            automargin: true
        },
    };

    Plotly.newPlot("ploidy_plot", data, layout);

};

// somatic variants plots

function build_variant_classification_plot(){

    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                // keep entries pertaining to vc_list
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                df_list.push(wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)})['Variant_Classification']);
            }
        }
    }
    // variant classification count
    let counts = dfd.concat({ df_list: df_list, axis: 0 }).value_counts()

    let trace1 = {
        x: counts['index_arr'],
        y: counts['data'],
        type: 'bar',
        transforms: [{
            type: 'sort',
            target: 'y',
            order: 'descending'
          }]
    };

    let data = [trace1];

    let layout = {
        title: 'Variant Classificaiton',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Variant Classification' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("variant_class_plot", data, layout);

}

function build_variant_type_plot(){

    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                // keep entries pertaining to vc_list
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                df_list.push(wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)})['Variant_Type']);
            }
        }
    }
    // variant type count
    let counts = dfd.concat({ df_list: df_list, axis: 0 }).value_counts()

    let trace1 = {
        x: counts['index_arr'],
        y: counts['data'],
        type: 'bar',
        transforms: [{
            type: 'sort',
            target: 'y',
            order: 'descending'
          }]
    };

    let data = [trace1];

    let layout = {
        title: 'Variant Type',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Variant Type' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("variant_type_plot", data, layout);

}

function build_snv_class_plot(){

    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                let ref = wes_data[i]['somatic']['maf']['Reference_Allele']['data'];
                let allele2 = wes_data[i]['somatic']['maf']['Tumor_Seq_Allele2']['data'];
                let class_list = [];
                // reference_allele > tumor_seq_allele2
                for (let j=0;j<ref.length; ++j) {
                    // bases only
                    if (checkAvailability(bases, ref[j]) && checkAvailability(bases, allele2[j])) {
                        // use snv_conv function to convert bases
                        class_list.push(snv_conv[ref[j]+'>'+allele2[j]])
                    }
                }
                df_list.push(new dfd.Series(class_list));
            }
        }
    }
    let counts = dfd.concat({ df_list: df_list, axis: 0 }).value_counts()

    let trace1 = {
        x: counts['index_arr'],
        y: counts['data'],
        type: 'bar',
        transforms: [{
            type: 'sort',
            target: 'y',
            order: 'descending'
          }]
    };

    let data = [trace1];

    let layout = {
        title: 'SNV Class',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'SNV Class' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("snv_class_plot", data, layout);

}


function build_variant_per_sample_plot(){

    let sample_ids = [];
    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], x: sample_ids, type: 'bar',name: vc_types[v]}));

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                vc = wes_data[i]['somatic']['maf']['Variant_Classification'].value_counts();
                // collect vc counts
                for (let j = 0; j < vc_types.length; ++j) {
                    data[j]['y'].push(vc['data'][vc['index_arr'].indexOf(vc_types[j])]);
                }
                // sample by sample
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }

    let layout = {
        title: 'Variants Per Sample',
        barmode: 'stack',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("variants_per_sample_plot", data, layout);

}


function build_variant_classification_summary_plot(){

    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], type: 'box',name: vc_types[v]}));

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                // vc count per sample
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                let vc_counts = wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)})['Variant_Classification'].value_counts();
                for (let j = 0; j < vc_types.length; ++j) {
                    data[j]['y'].push(vc_counts['data'][vc_counts['index_arr'].indexOf(vc_types[j])]);
                }
            }
        }
    }

    let layout = {
        title: 'Variant Classification Summary',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Variant Classification' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
        //legend: {traceorder: 'reversed'},
        showlegend: false
    };

    Plotly.newPlot("variant_class_summary_plot", data, layout);

}

function build_top_10_genes_plot(){
    // keep entries pertaining to vc_list
    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                df_list.push(wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)}));
            }
        }
    }
    // count top 10 gene occurances, all samples pooled
    let combined_df = dfd.concat({ df_list: df_list, axis: 0 })
    let gene_counts = combined_df['Hugo_Symbol'].value_counts()
    let keep_idx = findIndicesOfMax(gene_counts['data'], 20);
    let top_10_genes = keep_idx.map(i => gene_counts['index_arr'][i]);
    
    let combined_df_hs = combined_df['Hugo_Symbol']['data'];
    df_list = [];
    // keep entries pertaining to top_10_genes
    for (let i = 0; i < top_10_genes.length; ++i) {
        let keep_idx = [];
        for (let j=0;j<combined_df_hs.length; ++j) {
            if (combined_df_hs[j] == top_10_genes[i]) {
                keep_idx.push(j)
            }
        }
        df_list.push(combined_df.iloc({rows: Array.from(keep_idx)}));
    }

    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], x: top_10_genes, type: 'bar',name: vc_types[v]}));
    // vc counts by hugo_symbol
    for (let i = 0; i < top_10_genes.length; ++i) {
        vc = df_list[i]['Variant_Classification'].value_counts();
        for (let j = 0; j < vc_types.length; ++j) {
            data[j]['y'].push(vc['data'][vc['index_arr'].indexOf(vc_types[j])]);
        }
    }

    let layout = {
        title: 'Top 10 Mutated Genes',
        barmode: 'stack',
        xaxis: {
            title: { text: 'Gene' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("top_10_genes_plot", data, layout);

}