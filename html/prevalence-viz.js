/**
 * Generate the prevalence posterior pdf given the specified
 * parameters.
 * 
 * @param {number} n         - subjects tested
 * @param {number} k         - subjects that tested positive
 * @param {number} n_u       - known negative validation subjects
 * @param {number} k_u       - false positive validation subjects
 * @param {number} n_v       - known positive validation subjects
 * @param {number} k_v       - true  positive validation subjects
 * @param {number} [N=1000]  - monte carlo sample size 
 * @param {number} [M=1000]  - posterior grid size
 * @param {number} [a_u=1]   - prior parameter for false positive rate u: Beta(a_u,b_u)
 * @param {number} [b_u=1]   - prior parameter for false positive rate u: Beta(a_u,b_u)
 * @param {number} [a_v=1]   - prior parameter for true  positive rate v: Beta(a_v,b_v)
 * @param {number} [b_v=1]   - prior parameter for true  positive rate v: Beta(a_v,b_v)
 * @param {number} [w=1]     - population weighting
 * @returns {number[]} the posterior pdf
 */
function posterior({n, k, n_u, k_u, n_v, k_v, N=1000, M=1000, a_u=1, b_u=1, a_v=1, b_v=1, w=1}) {
    // reuse u,v samples for all posterior (theta) values
    let us=[], vs=[], duvs=[];
    for (let i=0; i<N; i++) { 
	let u=0, v=0;
	while (u >= v) {
            u = jStat.beta.sample(k_u + a_u, n_u - k_u + b_u);
            v = jStat.beta.sample(k_v + a_v, n_v - k_v + b_v);
	}
	us.push(u);
	vs.push(v);
	// to ameliorate numerical issues when Bu, Bv ~ 1, use
	// B(u;a,b) = 1-B(1-u;b,a). Then B(v;a,b) - B(u;a,b) =
	// B(1-u;b,a) - B(1-v;b,a) but the latter is the difference of
	// two numbers close to zero rather than 1.
	let Bu  = jStat.ibeta(1-u, n-k+1, k+1);
	let Bv  = jStat.ibeta(1-v, n-k+1, k+1);
	let duv = (v-u)/(Bu - Bv)
	duvs.push(duv)
    }

    // theta values
    let pdf         = [];
    let total       = 0;
    let ct          = 0;
    for (let j=0; j<=M; j++) {
	let theta   = j/M;
	// serial for now
	let p_theta = 0;
	let cp      = 0;
	for (let i=0; i<N; i++) {
	    let u   = us[i];
	    let v   = vs[i];
	    let duv = duvs[i];
	    let p   = u + theta * (v-u);
	    let f   = jStat.beta.pdf(p, k+1, n-k+1);
	    // p_theta is compensated sum of f*duv
            let y   = f*duv - cp;
            let t   = p_theta + y; 
            cp      = (t - p_theta) - y;
            p_theta = t;
	}
	pdf.push(p_theta);

	// total is compensated sum of p_theta
        let y   = p_theta - ct;
        let t   = total + y; 
        ct      = (t - total) - y;
        total   = t;
    }

    // normalization
    total /= pdf.length
    for (let j=0; j<=M; j++) {
	pdf[j] /= total;
    }
    
    return pdf;
}

function parseParams() {
    let params = {};
    allParams.forEach((k) => params[k] = parseFloat(document.getElementById(k).value));
    return params;
}

function setFormParams(params) {
    allParams.forEach((k) => document.getElementById(k).value = params[k]);
}

// the chart
let chart = false;

// chart config
Chart.defaults.global.defaultFontSize = 12;
Chart.defaults.global.defaultFontColor = "#222";
Chart.defaults.global.defaultFontFamily = "roboto";
const default_config = {
    type: "line",
    data: {
	labels:   [],
	datasets: []
    },
    options: {
	responsive: true,
	aspectRatio: 1,
	title: {
	    display: true,
	    text: "Posterior Prevalence Probability",
	    fontSize: 18
	},
	tooltips: {
	    mode: "index",
	    intersect: false,
	    // show 2 decimal places
            callbacks: {
                label: function(tooltipItem, data) {
                    let label = data.datasets[tooltipItem.datasetIndex].label || "";
                    if (label) {
                        label += ": ";
                    }
                    label += tooltipItem.yLabel.toFixed(2);
                    return label;
                }
            }
	},
	hover: {
	    mode: "nearest",
	    intersect: true
	},
	legend: {
	    display: false
	},
	annotation: {
            annotations: []
        },
	scales: {
	    xAxes: [{
		display: true,
		scaleLabel: {
		    display: true,
		    labelString: "prevalence(%)",
		    fontSize: 14,
		    fontStyle: "bold"
		}
	    }],
	    yAxes: [{
		display: true,
		scaleLabel: {
		    display: true,
		    labelString: "density",
		    fontSize: 14
		}
	    }]
	}
    }
};

function summaryStats(data, alpha) {
    // credible interval, median, and mode
    let lower=-1, upper=-1, median=-1, mode=-1, max=-1;
    let len = data.length;
    let cumsum = data.reduce((a, x, i) => [...a, x/len + (a[i-1] || 0)], []);
    cumsum.unshift(0);
    for (let i=0; i<cumsum.length; i++) {
	if (lower == -1 && cumsum[i] > alpha/2) {
	    lower = i;
	    // interpolate
	    if (i > 0) lower -= (cumsum[i] - alpha/2) / (cumsum[i] - cumsum[i-1]);
	}
	if (upper == -1 && cumsum[i] > 1 - alpha/2) {
	    upper = i;
	    // interpolate
	    if (i > 0) upper -= (cumsum[i] - (1-alpha/2)) / (cumsum[i] - cumsum[i-1]);
	}
	if (median == -1 && cumsum[i] > 1/2) {
	    median = i;
	    // interpolate
	    if (i > 0) median -= (cumsum[i] - 1/2) / (cumsum[i] - cumsum[i-1]);
	}
	if (data[i] > max) {
	    max = data[i];
	    mode = i;
	}
    }
    return [lower, upper, median, mode];
}

// annotations
const annotation = {
    drawTime: "afterDatasetsDraw",
    id: "",           // overridden
    type: "line",
    mode: "vertical",
    scaleID: "x-axis-0",
    value: -1,        // overridden
    borderColor: "black",
    borderWidth: 2,
    borderDash: [2, 2],
    label: {
	backgroundColor: "grey",
	content: "",  // overridden
	yAdjust: 0,   // overridden
	xAdjust: 0,   // overridden
	rotation: 90, // broken
	enabled: true
    }
};

function makeAnnotation(name, value, offset, len, weight, xAdjust, yAdjust) {
    let ann           = Object.assign({}, annotation);
    ann.label         = Object.assign({}, annotation.label);
    ann.id            = name;
    // position on the x-axis
    ann.value         = value - offset;
    // label corresponding to the position on the x-axis
    let mapped        = value * 100 * weight / len;
    ann.label.content = [`${name}:${mapped.toFixed(2)}%`];
    ann.label.xAdjust = xAdjust;
    ann.label.yAdjust = yAdjust;
    return ann;
}

// citation annotation - fake line with label
let citation = {
    drawTime: "afterDatasetsDraw",
    id: "citation",
    type: "line",
    mode: "horizontal",
    scaleID: "y-axis-0",
    value: 0, // overridden
    borderColor: "transparent",
    borderWidth: 0,
    label: {
	backgroundColor: "transparent",
	content: "Citation: Baxter, Bayesian beta-binomial prevalence estimation using an imperfect test (2020)",
	xAdjust: -40,
	fontSize: 10,
	fontColor:"#AAA",// "#FF4500",
	enabled: true
    }    
};

function updateUrl() {
    let paramFormElem = document.getElementById(paramForm);
    let search = new URLSearchParams(new FormData(paramFormElem));
    let url = new URL(window.location.href);
    url.search = search.toString();
    let params = parseParams();
    history.pushState(params, "Posterior Prevalence Probability (PPP)", url);
}

const updateChart = async () => {
    // read the parameters
    let params = parseParams();
    
    // calculate the prevalence posterior
    let data = posterior(params);

    let alpha = 1 - params.CI / 100;

    // credible interval, median, and mode
    let [lower, upper, median, mode] = summaryStats(data, alpha);
    
    // discard neglible pdf regions adjacent to boundary
    let max = Math.max(...data);
    let thresh = max / 1000;
    let len = data.length;
    let begin = 0, end = len-1;
    while (begin < len && begin < lower && data[begin] < thresh) ++begin;
    while (end > begin && end > upper   && data[end]   < thresh) --end;
    data = data.slice(begin, end+1);

    // retrieve or create the chart config
    let config = chart ? chart : default_config;
    
    // replace labels and data
    let weight  = params.w;
    config.data.labels.length = 0;
    for (let i=begin; i<=end; i++) {
	let x = 100 * i * weight / len;
	if (x > 100) x = 100;
	config.data.labels.push(x.toFixed(1));
    }
    let dataset = {
	backgroundColor: "#2a52be"/*"#69b3a2"*/,
	borderColor: "#000",
	borderWidth: 2,
	data: data,
	pointRadius:0,
	fill: true
    };
    let first = (config.data.datasets.length == 0);
    config.data.datasets.length=0;
    config.data.datasets.push(dataset);

    // annotations
    if (chart) config.annotation.elements = []; // hack to clear computed annotations
    let annotations = [];
    let yAdjust = 25
    annotations.push(makeAnnotation(`lower ${(100*alpha/2).toFixed(1)}%`     , lower , begin, len, weight, 55, 0        ));
    annotations.push(makeAnnotation(`upper ${(100*(1-alpha/2)).toFixed(1)}%` , upper , begin, len, weight, 55, yAdjust  ));
    annotations.push(makeAnnotation("median"                                 , median, begin, len, weight, 45, 2*yAdjust));
    annotations.push(makeAnnotation("mode"                                   , mode  , begin, len, weight, 40, 3*yAdjust));
    citation.value = max/20;
    annotations.push(citation);
    config.options.annotation.annotations = annotations;
    
    // update
    if (first) {
	// create the chart
	let ctx = document.getElementById("prevalence_viz").getContext("2d");
	chart = new Chart(ctx, config);
    } else {
	// update the chart
	chart.update();
    }
}

// event listeners
document.getElementById("paramForm").addEventListener(
    "submit",
    (event) => {
	event.preventDefault();
	updateUrl();
	updateChart();
    });

window.addEventListener(
    "popstate",
    (event) => {
	setFormParams(history.state);
	updateChart();
    });
