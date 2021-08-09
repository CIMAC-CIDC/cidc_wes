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
    build_ti_tv_plot();
    build_lego_plot();
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

var ti_tv_conv = {
    "T>C":"Ti",
    "C>T":"Ti",
    "T>A":"Tv",
    "T>G":"Tv",
    "C>A":"Tv",
    "C>G":"Tv"
}
  
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

function sumObjectsByKey(...objs) {
    return objs.reduce((a, b) => {
      for (let k in b) {
        if (b.hasOwnProperty(k))
          a[k] = (a[k] || 0) + b[k];
      }
      return a;
    }, {});
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

    Plotly.newPlot("GC_Content_plot", data, layout);

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
    };

    let data = [trace1];

    let layout = {
        title: 'SNV Class',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'SNV Class' },
            automargin: true,
            categoryorder: "array",
            categoryarray:  ['T>G','T>A','T>C','C>T','C>G','C>A']
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("SNV_Class_plot", data, layout);

}


function build_variant_per_sample_plot(){

    let sample_ids = [];
    let total_counts = [];
    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], x: sample_ids, type: 'bar',name: vc_types[v]}));

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) { 
                vc = wes_data[i]['somatic']['maf']['Variant_Classification'].value_counts();
                let total_count = 0
                // collect vc counts
                for (let j = 0; j < vc_types.length; ++j) {
                    data[j]['y'].push(vc['data'][vc['index_arr'].indexOf(vc_types[j])]);
                    if (typeof vc['data'][vc['index_arr'].indexOf(vc_types[j])] != 'undefined'){
                        total_count = total_count + vc['data'][vc['index_arr'].indexOf(vc_types[j])];
                    }
                }
                // sample by sample
                sample_ids.push(wes_data[i]['id']);
                total_counts.push(total_count);
            }
        }
    }
    // sort sample ids by total counts
    let indices = [...total_counts.keys()].sort((a, b) => total_counts[b] - total_counts[a]);
    let sample_ids_reordered = [sample_ids].map(a => indices.map(i => a[i]))[0];

    let layout = {
        title: 'Variants Per Sample',
        barmode: 'stack',
        xaxis: {
            title: { text: 'Run' },
            automargin: true,
            categoryorder: 'array',
            categoryarray: sample_ids_reordered
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

function build_ti_tv_plot(){

    let snv = {'T>G':[],'T>A':[],'T>C':[],'C>T':[],'C>G':[],'C>A':[]}
    let df_list = []
    let stack_colors = ['#ff7f0e','#d62728','#8c564b','#7f7f7f','#17becf','#ff7f0e'];
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

    let ti_tv = {'Ti':Array(df_list.length).fill(0),'Tv':Array(df_list.length).fill(0)}

    for (let i = 0; i < df_list.length; ++i) {
        let counts = df_list[i].value_counts()
        let counts_sum = counts['data'].reduce((a, b) => a + b, 0)
        for (let j = 0; j < Object.keys(snv).length; ++j) {
            if (checkAvailability(counts['index_arr'], Object.keys(snv)[j])) {
                snv[counts['index_arr'][j]].push((counts['data'][j]/counts_sum)*100);
                ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] = ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] + (counts['data'][j]/counts_sum)*100;
            } else {
                snv[counts['index_arr'][j]].push(0);
                ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] = ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] + 0;
            }
        }
    }
    let data = [];

    for (let i = 0; i < Object.keys(snv).length; ++i) {
        // snv box plot
        data.push(
            {
                name: Object.keys(snv)[i],
                y: Object.values(snv)[i],
                type: 'box',
                marker: { color: stack_colors[i] },
            }
        )
        // ti-tv boxplot
        data.push(
            {
                name: Object.keys(snv)[i],
                x: current_samples,
                y: Object.values(snv)[i],
                type: 'bar',
                xaxis: 'x3',
                yaxis: 'y3',
            }
        )
    }
// stacked barplot
    for (let i = 0; i < Object.keys(ti_tv).length; ++i) {
        data.push(
            {
                name: Object.keys(ti_tv)[i],
                y: Object.values(ti_tv)[i],
                type: 'box',
                xaxis: 'x2',
                yaxis: 'y2',
            }
        )
    }

    let layout = {
        title: 'Ti-Tv',
        barmode: 'stack',
        showlegend: false,
        // subplots
        xaxis: {domain: [0, 0.7], anchor: 'y'},
        yaxis: {domain: [0.55, 1], anchor: 'x', range: [0, 100], title: { text: '% Mutations' }},
        yaxis2: {domain: [0.55, 1], anchor: 'x2', range: [0, 100]},
        xaxis2: {domain: [0.8, 1], anchor: 'y2'},
        xaxis3: {domain: [0, 1], anchor: 'y3'},
        yaxis3: {domain: [0, 0.45], anchor: 'x3', range: [0, 100], title: { text: '% Mutations' }},
    };

    Plotly.newPlot("Ti-Tv_plot", data, layout);

}

function build_lego_plot(){
    // each tri
    let tri_template = {
        "A_A": 0,
        "A_C": 0,
        "A_G": 0,
        "A_T": 0,
        "C_A": 0,
        "C_C": 0,
        "C_G": 0,
        "C_T": 0,
        "G_A": 0,
        "G_C": 0,
        "G_G": 0,
        "G_T": 0,
        "T_A": 0,
        "T_C": 0,
        "T_G": 0,
        "T_T": 0
    }
    // each snv
    let tri_data = {
        "C>A": { ...tri_template },
        "C>G": { ...tri_template },
        "C>T": { ...tri_template },
        "T>A": { ...tri_template },
        "T>C": { ...tri_template },
        "T>G": { ...tri_template }
    }

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['tri_matrix'] != undefined) {
                for (var snv in tri_data) {
                    if (Object.prototype.hasOwnProperty.call(wes_data[i]['somatic']['tri_matrix'], snv)) {
                        // add to tri_data
                        tri_data[snv] = sumObjectsByKey({ ...tri_data[snv]},
                            { ...wes_data[i]['somatic']['tri_matrix'][snv]})
                    }
                }
            }
        }
    }
    // normalize tri_data in each snv
    for (var snv in tri_data) {
        let tri_sum = Object.values(tri_data[snv]).reduce((a, b) => a + b, 0)
            Object.keys(tri_data[snv])
                .forEach(function (tri) {
                    tri_data[snv][tri] = (tri_data[snv][tri] / tri_sum) * 100
                });
        }

    let data = [{
        x: Object.keys(tri_data['C>A']),
        y: Object.values(tri_data['C>A']),
        type: 'bar',
        xaxis: 'x',
        yaxis: 'y',
        name: 'C>A'
    }, {
        x: Object.keys(tri_data['C>G']),
        y: Object.values(tri_data['C>G']),
        type: 'bar',
        xaxis: 'x2',
        yaxis: 'y2',
        name: 'C>G'
    }, {
        x: Object.keys(tri_data['C>T']),
        y: Object.values(tri_data['C>T']),
        type: 'bar',
        xaxis: 'x3',
        yaxis: 'y3',
        name: 'C>T'
    }, {
        x: Object.keys(tri_data['T>A']),
        y: Object.values(tri_data['T>A']),
        type: 'bar',
        xaxis: 'x4',
        yaxis: 'y4',
        name: 'T>A'
    }, {
        x: Object.keys(tri_data['T>C']),
        y: Object.values(tri_data['T>C']),
        type: 'bar',
        xaxis: 'x5',
        yaxis: 'y5',
        name: 'T>C'
    }, {
        x: Object.keys(tri_data['T>G']),
        y: Object.values(tri_data['T>G']),
        type: 'bar',
        xaxis: 'x6',
        yaxis: 'y6',
        name: 'T>G'
    },
    ];

    let layout = {
        title: 'Lego Plot',
        height: 450,
        width: 1044,
        // subplots
        xaxis: {domain: [0, 0.3], anchor: 'y'},
        yaxis: {domain: [0.58, 1], anchor: 'x', title: { text: '% Mutations' }},
        yaxis2: {domain: [0.58, 1], anchor: 'x2'},
        xaxis2: {domain: [0.35, .65], anchor: 'y2'},
        xaxis3: {domain: [.7, 1], anchor: 'y3'},
        yaxis3: {domain: [.58, 1], anchor: 'x3'},
        xaxis4: {domain: [0, 0.3], anchor: 'y4'},
        yaxis4: {domain: [0, .42], anchor: 'x4', title: { text: '% Mutations' }},
        yaxis5: {domain: [0, .42], anchor: 'x5'},
        xaxis5: {domain: [0.35, .65], anchor: 'y5'},
        xaxis6: {domain: [.7, 1], anchor: 'y6'},
        yaxis6: {domain: [0, .42], anchor: 'x6'},
    };

    Plotly.newPlot("lego_plot_plot", data, layout);

}