//change to true for init options
let show_unclassified = false;
//default value for k_mer_threshold
let k_mer_threshold = 0.01;
//set after init when df is given
let species_list = [];
let checkbox_list = [];
// dataframe of thresholdXspecies
let df = "";

const species_colors = {};
const TOP_SPECIES = 10; // should be if df has more columns.
let actual_top_species = 0; // will be set on first run.
let non_contaminant_col_name = "";

const initResultsScript = (data, non_contaminant_col_name) => {
    const json_data = JSON.parse(data);
    non_contaminant_col_name = non_contaminant_col_name;
    runResultsScript(json_data)
};

//set functions
const update_species_list = (item) => {

    if (item.checked) {
        if (!species_list.includes(item.id)) {
            species_list.push(item.id)
        }
    } else {
        if (species_list.includes(item.id)) {
            item_index = species_list.indexOf(item.id);
            species_list.splice(item_index, 1);
        }
    }
    console.log(species_list)
    document.getElementById("species_list").value = [...species_list];
}
const update_k_mer_threshold = (new_val) => {
    k_mer_threshold = new_val;
    document.getElementById("k_mer_threshold").value = new_val;
}

const runResultsScript = (json_data) => {    
    let new_index = []
    let new_data = {}
    // const json_data = JSON.parse('{{ data|safe }}');
    Object.entries(json_data).map(([key, val]) => {
        new_data[key] = Object.values(val);
        new_index = Object.keys(val);
    })
    new_index = new_index.map((item) => parseFloat(item))
    df = new dfd.DataFrame(new_data, {index: new_index});
    
    //updates when pressing checkboxes
    df.columns.forEach(item => update_species_list({id: item, checked: true}))
    update_species_list({id:non_contaminant_col_name, checked: false});
    actual_top_species = Math.min(df.columns.length-1, TOP_SPECIES)
    species_list.forEach((item,idx) => {
        species_colors[item] = colors[idx]
    })

    update_k_mer_threshold(k_mer_threshold)    

    Chart.register("chartjs-plugin-annotation");
    bar_chart, pie_chart = createResultsCharts()
    
    draw_kmer_hist(bar_chart, pie_chart)
    draw_species_pie_chart(pie_chart)
    draw_species_list()
}

const createResultsCharts = () => {
    bar_chart = new Chart('bar_chart', {
        type: "bar",
        data: {},
        options: {}
    });

    pie_chart = new Chart('pie_chart', {
        type: "doughnut",
        data: {},
        options: {}
    });
    return bar_chart, pie_chart;
}

const draw_species_pie_chart = (chart) => {
    
    let sum_total = df.sum({axis: 1}).sum({axis: 0});
    const threshold_list = df.index.filter(item => item >= k_mer_threshold)

    let sum_classified_df = df.loc({ rows: threshold_list, columns: species_list }).sum({axis: 1});

    let sum_unclassified = sum_total - sum_classified_df.sum({axis: 0});
    let colors_classified = species_list.map((item) => species_colors[item]);

    const sorted_by_freq = get_sorted_species(sum_classified_df);
    const top_species = sorted_by_freq.slice(0, actual_top_species)
    const top_num_reads = sum_classified_df.loc(sorted_by_freq.slice(0,actual_top_species)).values
    const other_species_num_reads = sum_classified_df.loc(sorted_by_freq.slice(actual_top_species,)).sum()

    chart.data.labels = [...top_species, "other"]
    show_unclassified ? chart.data.labels.push("unclassified") : null;

    chart.data.datasets = [{
        label: 'Classified',
        data: [...top_num_reads, other_species_num_reads],
        hoverOffset: 4,
        backgroundColor: [...colors_classified, "orange"]
        }
    ];
    show_unclassified ? chart.data.datasets[0].data.push(sum_unclassified) : null;
    show_unclassified ? chart.data.datasets[0].backgroundColor.push("blue") : null;


    chart.options.borderRadius = 10
    chart.options.radius = "75%"
    chart.options.plugins.legend.display = false;

    chart.update()
}

const draw_kmer_hist = (chart, pie_chart) => {

    const sum_total_df = df.sum({axis: 0});
    const y_max = sum_total_df.max();
    let sum_classified_df = df.loc({ columns: species_list }).sum({axis: 0});
    const classified_indexes = sum_classified_df.index
    const classified_values = classified_indexes.map((index) => {
        return index < k_mer_threshold ? 0 : sum_classified_df.loc([index]).values[0];
    })
    sum_classified_df = new dfd.Series(classified_values, { index: classified_indexes });
    let sum_unclassified_df = sum_total_df.sub(sum_classified_df)
    const bins_x = sum_total_df.index
    const classified_y = sum_classified_df.values
    const unclassified_y = sum_unclassified_df.values


    chart.data.labels = bins_x
    chart.data.datasets = [
        {
            data: classified_y,
            backgroundColor: "orange",
            label: 'Classified'

        }
    ];
    
    show_unclassified ? chart.data.datasets.push({
            data: unclassified_y,
            backgroundColor: "blue",
            label: 'Unclassified'
        }) : null;


    chart.options = {
        onClick: (event) => { 

            const xTop = chart.chartArea.left;
            const xBottom = chart.chartArea.right;
            const xMin = chart.scales.x.min;
            const xMax = chart.scales.x.max;
            let newX = 0;

            if (event.x <= xBottom && event.x >= xTop) {
                newX = Math.abs((event.x - xTop) / (xBottom - xTop));
                newX = newX * (Math.abs(xMax - xMin)) + xMin;
            };
            update_k_mer_threshold(newX);
            draw_kmer_hist(chart, pie_chart)
            draw_species_pie_chart(pie_chart)
            draw_species_list()
        },
        plugins: {
            title: {
                display: true,
                text: 'classification by threshold k-mer',
                font: {
                    size: 24
                }
            },
            annotation: {
                annotations: {
                    line: {
                        type: 'line',
                        scaleID: 'x',
                        value: parseFloat(k_mer_threshold),
                        borderColor: 'rgb(255, 99, 132)',
                        borderWidth: 2,
                    }
                }
            },
            legend: {
                display: false
            }
        },
        scales: {
            x: {
                beginAtZero: true,
                stacked: true,
                type: 'linear',
                min: 0.0,
                max: 1.0,
                offset: false,
                title: {
                    display: true,
                    text: 'Percent Threshold',
                    font: {
                        size: 16
                    }
                }
            },
            y: {
                stacked: true,
                type: 'logarithmic',
                title: {
                    display: true,
                    text: 'Number of reads (log-scale)',
                    font: {
                        size: 16
                    }
                }
            }
        }
    };
    chart.update()
}

const draw_species_list = () => {
    const species_container = document.getElementById("species_container");
    const threshold_list = df.index.filter(item => item >= k_mer_threshold);
    let sum_classified_df = df.loc({ rows: threshold_list}).sum({axis: 1});
    
    const sorted_by_freq = get_sorted_species(sum_classified_df);
    
    const top_chosen_species = sorted_by_freq.slice(0,actual_top_species).map((item) => {
        let toggle = document.createElement("label");
        toggle.setAttribute("id", "toggle_" + item);
        const checkbox = document.createElement("input");
        checkbox.setAttribute("type", "checkbox");
        checkbox.setAttribute("id", item);
        checkbox.setAttribute("class", "h-5 w-5")
        checkbox.checked = species_list.includes(item) ? true : false;
        checkbox.addEventListener("click", (change) => {
            const item = change.target;
            update_species_list(item);
            //update the species_list to sent to backend
            draw_kmer_hist(bar_chart, pie_chart);
            draw_species_pie_chart(pie_chart);
        });

        toggle.appendChild(checkbox);
        let label = document.createElement("span");
        label.innerText = item;
        toggle.appendChild(label);
        return toggle;
    });

    const other_species = sorted_by_freq.slice(actual_top_species,).map((item) => {
        const listitem = document.createElement("li");
        const toggle = document.createElement("label");
        listitem.setAttribute("id", "toggle_" + item);
        const checkbox = document.createElement("input");
        checkbox.setAttribute("type", "checkbox");
        checkbox.setAttribute("id", item);
        checkbox.setAttribute("class", "h-5 w-5 mx-5")
        checkbox.checked = species_list.includes(item) ? true : false;

        checkbox.addEventListener("click", (change) => {
            const item = change.target;
            update_species_list(item);
            draw_kmer_hist(bar_chart, pie_chart);
            draw_species_pie_chart(pie_chart);
        });
        let label = document.createElement("span");
        label.innerText = item;
        toggle.appendChild(checkbox)
        toggle.appendChild(label);
        listitem.appendChild(toggle);
        return listitem;
    });

    const other_species_container = document.createElement("ul");
    const listitem = document.createElement("li");
    
    const label = document.createElement("span");
    label.innerText = "Other Organisms";
    const other_species_toggle = document.createElement("input");
    const other_species_label = document.createElement("label");

    other_species_toggle.setAttribute("type", "checkbox");
    other_species_toggle.setAttribute("value", "other");
    other_species_toggle.setAttribute("id", "toggle_other_species");
    other_species_toggle.setAttribute("class", "h-5 w-5");
    other_species_toggle.checked = true;

    other_species_toggle.addEventListener("click", (change) => {
        sorted_by_freq.slice(actual_top_species,).forEach((item) => {

            document.getElementById(item).checked = change.target.checked
            update_species_list({id: item, checked: change.target.checked})

        });
        draw_kmer_hist(bar_chart, pie_chart);
        draw_species_pie_chart(pie_chart);
    });

    other_species_label.appendChild(other_species_toggle)
    other_species_label.appendChild(label)


    listitem.appendChild(other_species_label);
    const other_species_list = document.createElement("ul");
    other_species_list.replaceChildren(...other_species)
    listitem.appendChild(other_species_list)
    other_species_container.replaceChildren(listitem);

    species_container.children[0].replaceChildren(...top_chosen_species)
    species_container.children[1].replaceChildren(other_species_container)

}

const toggle_unclassified = () => {
    show_unclassified = !show_unclassified;
    draw_kmer_hist(bar_chart, pie_chart)
    draw_species_pie_chart(pie_chart)
}

// sorted species list by reads without unclassified entry
const get_sorted_species = (df) => {
    let sort_index = df.sort_values({"ascending": false }).index
    let sorted_by_freq = df.iloc(sort_index).index
    sorted_by_freq = sorted_by_freq.filter(item => item != non_contaminant_col_name)
    return sorted_by_freq;
}