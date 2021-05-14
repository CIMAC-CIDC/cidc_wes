function checkAvailability(arr, val) {
    return arr.some(function (arrVal) {
        return val === arrVal;
    });
}

function build_plot(){
    build_mapping_plot();
    build_coverage_plot();
    build_mean_quality_plot();
    build_gc_content_plot();
    build_insert_size_plot();
    build_clonality_plot();
    build_purity_plot();
    build_ploidy_plot();
}

// Data quality plots

var sample_type = ['normal', 'tumor'];
var sample_type_suffix = ['.N', '.T'];

function build_mapping_plot() {

    var total_reads = [];
    var mapped_reads = [];
    var dedup_reads = [];
    var sample_ids = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    total_reads.push(wes_data[i][sample_type[j]]['alignment']['total_reads'])
                    mapped_reads.push(wes_data[i][sample_type[j]]['alignment']['mapped_reads'])
                    dedup_reads.push(wes_data[i][sample_type[j]]['alignment']['dedup_reads'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }

    var trace1 = {
        x: sample_ids,
        y: total_reads,
        name: 'Total Reads',
        type: 'bar'
    };

    var trace2 = {
        x: sample_ids,
        y: mapped_reads,
        name: 'Mapped Reads',
        type: 'bar'
    };

    var trace3 = {
        x: sample_ids,
        y: dedup_reads,
        name: 'Dedup Reads',
        type: 'bar'
    };

    var data = [trace1, trace2, trace3];

    var layout = {
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

    var mean_depth = [];
    var sample_ids = [];

    var hover_text = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    mean_depth.push(wes_data[i][sample_type[j]]['coverage']['mean_depth'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
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

    var trace1 = {
        y: sample_ids,
        x: mean_depth,
        name: 'Mean Depth',
        type: 'bar',
        orientation: 'h',
        hovertemplate: '%{text}<extra></extra>',
        text: hover_text,
    };

    var data = [trace1];

    var layout = {
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

    var gc_content = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    gc_content.push({
                        // x: [...Array(101).keys()],
                        y: wes_data[i][sample_type[j]]['alignment']['gc_content'],
                        mode: 'lines',
                        type: 'scatter',
                        name: wes_data[i]['id'] + sample_type_suffix[j]
                    })
                }
            }
        }
    }

    var data = gc_content;

    var layout = {
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

    var insert_size = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    insert_size.push({
                        // x: [...Array(101).keys()],
                        y: wes_data[i][sample_type[j]]['alignment']['insert_size'],
                        mode: 'lines',
                        type: 'scatter',
                        name: wes_data[i]['id'] + sample_type_suffix[j]
                    })
                }
            }
        }
    }

    var data = insert_size;

    var layout = {
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

    var mean_quality = [];
    var sample_ids = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    mean_quality.push(wes_data[i][sample_type[j]]['alignment']['mean_quality_score'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }

    var trace1 = {
        x: sample_ids,
        y: mean_quality,
        type: 'bar'
    };

    var data = [trace1];

    var layout = {
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

    var clonality = [];
    var sample_ids = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                clonality.push(wes_data[i]['copy_number']['clonality'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    var trace1 = {
        x: sample_ids,
        y: clonality,
        type: 'bar'
    };

    var data = [trace1];

    var layout = {
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

    var purity = [];
    var sample_ids = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                purity.push(wes_data[i]['copy_number']['purity'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    var trace1 = {
        x: sample_ids,
        y: purity,
        type: 'bar'
    };

    var data = [trace1];

    var layout = {
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

    var ploidy = [];
    var sample_ids = [];

    for (i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                ploidy.push(wes_data[i]['copy_number']['ploidy'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    var trace1 = {
        x: sample_ids,
        y: ploidy,
        mode: 'lines',
        type: 'scatter',
    };

    var data = [trace1];

    var layout = {
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