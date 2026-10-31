// coi-sw.js — cross-origin isolation for a host that cannot send headers (GitHub Pages): SharedArrayBuffer, and with
// it the solver threads of js/parallel.js, exists only in an isolated page. Loaded by index.html as a plain script,
// first thing in <head>, this file registers itself as a service worker for the app's folder; the service worker adds
// Cross-Origin-Opener-Policy and Cross-Origin-Embedder-Policy to every response, and the page reloads once so that it
// runs isolated. Active only on *.github.io, or with ?coi=1 in the address (a local test; ?coi=0 removes it again):
// a host that sends the headers itself (multi_threading/serve.py, Cloudflare Pages with a _headers file) needs none,
// and localhost stays free of a service worker that would also catch other apps served on the same port.
// Without isolation the app still works: one thread, the WebAssembly kernels, the same results.
if (typeof window === 'undefined') {
  // ---------------------------------------------------------------- the service worker
  self.addEventListener('install', () => self.skipWaiting());
  self.addEventListener('activate', (e) => e.waitUntil(self.clients.claim()));
  self.addEventListener('fetch', (e) => {
    const req = e.request;
    if (req.cache === 'only-if-cached' && req.mode !== 'same-origin') return;      // leave those to the browser
    e.respondWith(fetch(req).then((res) => {
      if (res.status === 0) return res;                                            // opaque (cross-origin no-cors): as it is
      const h = new Headers(res.headers);
      h.set('Cross-Origin-Opener-Policy', 'same-origin');
      h.set('Cross-Origin-Embedder-Policy', 'require-corp');
      if (!h.has('Cross-Origin-Resource-Policy')) h.set('Cross-Origin-Resource-Policy', 'cross-origin');
      const empty = res.status === 101 || res.status === 204 || res.status === 205 || res.status === 304;
      return new Response(empty ? null : res.body, { status: res.status, statusText: res.statusText, headers: h });
    }));
  });
} else (() => {
  // ---------------------------------------------------------------- the page: register, reload once
  const me = document.currentScript && document.currentScript.src, sw = navigator.serviceWorker;
  if (!me || !sw || !window.isSecureContext) return;
  const want = new URLSearchParams(location.search).get('coi');
  if (want === '0') {                                                              // remove the service worker, load without it
    sw.getRegistrations().then((rs) => Promise.all(rs.filter((r) => [r.active, r.waiting, r.installing].some((w) => w && w.scriptURL === me)).map((r) => r.unregister())))
      .then((gone) => { if (gone.length) location.replace(location.pathname); });
    return;
  }
  let ss = null; try { ss = window.sessionStorage; } catch (e) { /* storage blocked: no reload, so no isolation */ }
  const KEY = 'coi-sw-reloaded';
  if (window.crossOriginIsolated) { if (ss) ss.removeItem(KEY); return; }
  if (!(location.hostname.endsWith('.github.io') || want === '1') || !ss || ss.getItem(KEY)) return;   // reloaded once already: one thread it is
  sw.register(me).then((reg) => {
    const reload = () => { ss.setItem(KEY, '1'); location.reload(); };
    if (sw.controller || (reg.active && !reg.installing && !reg.waiting)) reload();   // active, but this load came without the headers
    else sw.addEventListener('controllerchange', reload, { once: true });
  }, (err) => console.warn('coi-sw: no service worker:', err));
})();
