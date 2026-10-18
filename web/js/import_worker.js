// Worker for the mesh import (import_pipeline.js): read, repair, remesh and intrinsic Delaunay off the main thread.
// Message in: {id, buf (transferred), name, opts}. Out: {id, progress, text} while it runs, then {id, done, data}
// (the carrier arrays transferred) or {id, error}.
import { prepareCarrier } from './import_pipeline.js';

self.onmessage = (e) => {
  const { id, buf, name, opts } = e.data;
  try {
    const d = prepareCarrier(buf, name, opts || {}, (progress, text) => self.postMessage({ id, progress, text }));
    const { stats, ...data } = d; void stats;
    const xfer = [data.V, data.F, data.Fi, data.Li, data.K, data.B].map(a => a.buffer);
    self.postMessage({ id, done: true, data }, xfer);
  } catch (err) { self.postMessage({ id, error: String(err && err.message || err) }); }
};
