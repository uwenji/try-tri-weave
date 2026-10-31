// Small sparse toolkit: CSR matrices, SpMV, Jacobi-preconditioned conjugate gradients, a binary heap and
// multi-source Dijkstra. Float64 on the CPU; the GPU solver uploads f32 copies.

export function csrFromCOO(n, rows, cols, vals) {
  // sums duplicates; returns {n, indptr: Int32Array, indices: Int32Array, data: Float64Array}
  const nnz = rows.length;
  const count = new Int32Array(n + 1);
  for (let k = 0; k < nnz; k++) count[rows[k] + 1]++;
  for (let i = 0; i < n; i++) count[i + 1] += count[i];
  const pos = count.slice(0, n);
  const ci = new Int32Array(nnz), cv = new Float64Array(nnz);
  for (let k = 0; k < nnz; k++) { const p = pos[rows[k]]++; ci[p] = cols[k]; cv[p] = vals[k]; }
  // sort each row by column and merge duplicates (stable: duplicates are summed in input order). Short rows by insertion
  // sort in place, long ones by a stable index sort — the same order, so the same sums, as sorting every row
  const indptr = new Int32Array(n + 1), indices = new Int32Array(nnz), data = new Float64Array(nnz);
  let t = 0;
  for (let i = 0; i < n; i++) {
    const s = count[i], e = count[i + 1];
    if (e - s <= 64) {
      for (let a = s + 1; a < e; a++) {
        const c = ci[a], v = cv[a]; let b = a - 1;
        while (b >= s && ci[b] > c) { ci[b + 1] = ci[b]; cv[b + 1] = cv[b]; b--; }
        ci[b + 1] = c; cv[b + 1] = v;
      }
    } else {
      const idx = Array.from({ length: e - s }, (_, k) => s + k).sort((a, b) => ci[a] - ci[b]);
      const c2 = idx.map(k => ci[k]), v2 = idx.map(k => cv[k]);
      for (let k = 0; k < idx.length; k++) { ci[s + k] = c2[k]; cv[s + k] = v2[k]; }
    }
    let last = -1;
    for (let k = s; k < e; k++) {
      if (ci[k] === last) data[t - 1] += cv[k];
      else { indices[t] = ci[k]; data[t] = cv[k]; t++; last = ci[k]; }
    }
    indptr[i + 1] = t;
  }
  return { n, indptr, indices: indices.slice(0, t), data: data.slice(0, t) };
}

export function spmv(A, x, y) {
  const { n, indptr, indices, data } = A;
  for (let i = 0; i < n; i++) {
    let s = 0.0;
    for (let k = indptr[i]; k < indptr[i + 1]; k++) s += data[k] * x[indices[k]];
    y[i] = s;
  }
  return y;
}

export function csrDiagonal(A) {
  const d = new Float64Array(A.n);
  for (let i = 0; i < A.n; i++) for (let k = A.indptr[i]; k < A.indptr[i + 1]; k++) if (A.indices[k] === i) d[i] += A.data[k];
  return d;
}

export function dot(a, b) { let s = 0.0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

/** Jacobi-PCG for SPD operators given as a function applyA(x, out). Returns {x, iters, relres}. applyA may return
 *  dot(x, out) summed in index order (it then equals dot(x, out) bit for bit), which saves one pass per iteration. */
export function pcg(applyA, diag, b, x0, { maxIter = 500, tol = 1e-8, projectMean = null } = {}) {
  const n = b.length;
  const x = x0 ? Float64Array.from(x0) : new Float64Array(n);
  const r = new Float64Array(n), z = new Float64Array(n), p = new Float64Array(n), Ap = new Float64Array(n);
  applyA(x, Ap);
  for (let i = 0; i < n; i++) r[i] = b[i] - Ap[i];
  if (projectMean) projectMean(r);
  const bnorm = Math.sqrt(dot(b, b)) || 1.0;
  for (let i = 0; i < n; i++) { z[i] = r[i] / diag[i]; p[i] = z[i]; }
  let rz = dot(r, z), it = 0, rn = Math.sqrt(dot(r, r));
  while (it < maxIter && rn > tol * bnorm) {
    const ret = applyA(p, Ap);                                              // an operator may return dot(p, Ap) itself
    const pAp = typeof ret === 'number' ? ret : dot(p, Ap);
    if (!(pAp > 0)) break;
    const alpha = rz / pAp;
    let rzNew = 0, rr = 0;
    if (projectMean) {
      for (let i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
      projectMean(r);
      for (let i = 0; i < n; i++) z[i] = r[i] / diag[i];
      rzNew = dot(r, z); rr = dot(r, r);
    } else {
      // one pass for the update, the preconditioner and both dot products (same operations in the same order as the
      // separate loops, so the same numbers)
      for (let i = 0; i < n; i++) {
        x[i] += alpha * p[i];
        const ri = r[i] - alpha * Ap[i], zi = ri / diag[i];
        r[i] = ri; z[i] = zi; rzNew += ri * zi; rr += ri * ri;
      }
    }
    const beta = rzNew / rz; rz = rzNew;
    for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    rn = Math.sqrt(rr); it++;
  }
  return { x, iters: it, relres: rn / bnorm };
}

// ------------------------------------------------------------------ binary heap (min) of (key, id)
export class MinHeap {
  constructor(cap = 1024) { this.keys = new Float64Array(cap); this.ids = new Int32Array(cap); this.size = 0; }
  push(key, id) {
    if (this.size === this.keys.length) { const k2 = new Float64Array(this.size * 2); k2.set(this.keys); this.keys = k2; const i2 = new Int32Array(this.size * 2); i2.set(this.ids); this.ids = i2; }
    let i = this.size++;
    while (i > 0) { const p = (i - 1) >> 1; if (this.keys[p] <= key) break; this.keys[i] = this.keys[p]; this.ids[i] = this.ids[p]; i = p; }
    this.keys[i] = key; this.ids[i] = id;
  }
  pop() {
    const key = this.keys[0], id = this.ids[0];
    const lastK = this.keys[--this.size], lastI = this.ids[this.size];
    let i = 0;
    for (;;) {
      let c = 2 * i + 1; if (c >= this.size) break;
      if (c + 1 < this.size && this.keys[c + 1] < this.keys[c]) c++;
      if (this.keys[c] >= lastK) break;
      this.keys[i] = this.keys[c]; this.ids[i] = this.ids[c]; i = c;
    }
    if (this.size > 0) { this.keys[i] = lastK; this.ids[i] = lastI; }
    return [key, id];
  }
}

/** Multi-source Dijkstra on a CSR graph with edge weights in A.data. Returns {dist, source, pred}. */
export function dijkstra(G, sources) {
  const n = G.n, dist = new Float64Array(n).fill(Infinity), src = new Int32Array(n).fill(-1), pred = new Int32Array(n).fill(-1);
  const heap = new MinHeap(Math.max(64, n));
  for (const s of sources) { dist[s] = 0; src[s] = s; heap.push(0, s); }
  while (heap.size) {
    const [d, i] = heap.pop();
    if (d > dist[i]) continue;
    for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) {
      const j = G.indices[k], nd = d + G.data[k];
      if (nd < dist[j]) { dist[j] = nd; src[j] = src[i]; pred[j] = i; heap.push(nd, j); }
    }
  }
  return { dist, source: src, pred };
}
