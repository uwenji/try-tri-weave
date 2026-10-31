// A helper thread of the multi-threaded CG (parallel.js): it gets the shared control block once (and answers 'ready'),
// then one message per solve, and computes the chunks of rows and vector entries the coordinator (global_worker.js)
// hands out. Runs as a browser module worker and as a Node worker_threads worker (test_parallel.mjs).
import { helperMain } from './parallel.js';

if (typeof self !== 'undefined' && typeof self.addEventListener === 'function') {
  helperMain((handle) => { self.onmessage = (e) => handle(e.data); }, (msg) => self.postMessage(msg));
} else {
  const { parentPort } = await import('node:worker_threads');
  helperMain((handle) => parentPort.on('message', handle), (msg) => parentPort.postMessage(msg));
}
