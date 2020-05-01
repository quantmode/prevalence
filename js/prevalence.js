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
 * @returns {number[]} the posterior pdf
 */
export function posterior({n, k, n_u, k_u, n_v, k_v, N=1000, M=1000, a_u=1, b_u=1, a_v=1, b_v=1}) {
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

    // normalize
    total /= pdf.length
    for (let j=0; j<=M; j++) {
	pdf[j] /= total;
    }
    
    return pdf;
}
