"""Serve the web app on this machine and to other devices on the same network.

    python serve.py                 # http  on port 8080, prints the URLs to open from this Mac and from other devices
    python serve.py --https         # https on port 8443 with a self-signed certificate (made once with openssl)
    python serve.py --port 9000     # a busy default port (another server still running) moves on to the next free one
    python serve.py --verbose       # log every request and every dropped connection, with the client's address

Why --https: browsers expose WebGPU only in a "secure context". http://localhost counts as secure, but
http://192.168.x.x does not, so a phone or iPad on the Wi-Fi gets the (slower) CPU solver over plain http.
Over https it gets WebGPU after you accept the certificate warning once (Safari: "Show details → visit this
website"; Chrome: "Advanced → proceed"). The published Claude artifact is the no-setup alternative: it is https.
macOS may ask once to allow incoming connections for Python — allow it, or other devices cannot connect; the
network check printed at start-up says what to switch on when they still cannot.
"""
from __future__ import annotations

import argparse
import errno
import http.server
import os
import socket
import ssl
import subprocess
import sys
import threading
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
APP = 'TryTriWeave'
PEM = 'localsh-selfsigned.pem'


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map,
                      '.js': 'text/javascript', '.mjs': 'text/javascript', '.lshc': 'application/octet-stream', '.wasm': 'application/wasm', '.json': 'application/json'}

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store')              # always refetch edited modules / re-exported carriers (no-cache still let Chrome reuse modules)
        super().end_headers()

    def log_message(self, fmt, *args):                              # quieter log: skip the per-module chatter
        line = fmt % args
        if not VERBOSE and ('"GET /js/' in line or '"GET /vendor/' in line or '"GET /carriers/' in line or 'favicon' in line):
            return
        super().log_message(fmt, *args)


VERBOSE = False


# a browser drops connections all the time (preconnects it never uses, reloads, the certificate warning, a tab closed
# mid-download of a carrier); macOS then reports EBADF / ECONNRESET / EPIPE in the handler thread. Nothing is wrong: stay quiet.
QUIET_ERRNO = {errno.EBADF, errno.ECONNRESET, errno.EPIPE, errno.ECONNABORTED, errno.ETIMEDOUT, errno.ENOTCONN, errno.ESHUTDOWN}


class Server(http.server.ThreadingHTTPServer):
    daemon_threads = True
    tls = None                                                      # ssl.SSLContext when serving https

    def finish_request(self, request, client_address):
        """https: look at the first byte. A TLS hello starts with 0x16; anything else is a browser that sent plain http
        to the https port (Safari does that for a typed address without https://) - send it to the https address."""
        if self.tls is None:
            return super().finish_request(request, client_address)
        request.settimeout(15)
        first = request.recv(1, socket.MSG_PEEK)
        if not first:
            return
        if first[0] != 0x16:
            head = request.recv(8192).decode('latin-1', 'replace')
            lines = head.split('\r\n')
            path = lines[0].split(' ')[1] if len(lines[0].split(' ')) > 1 else '/'
            host = next((l.split(':', 1)[1].strip() for l in lines[1:] if l.lower().startswith('host:')), '')
            if not host:
                host = f'{(lan_addresses() or ["localhost"])[0]}:{self.server_address[1]}'
            location = f'https://{host}{path}'
            body = f'<a href="{location}">{location}</a>\n'.encode()
            request.sendall(f'HTTP/1.1 301 Moved Permanently\r\nLocation: {location}\r\nContent-Type: text/html\r\n'
                            f'Content-Length: {len(body)}\r\nConnection: close\r\n\r\n'.encode() + body)
            if VERBOSE:
                print(f'{client_address[0]} - plain http on the https port: sent to {location}')
            return
        tls = self.tls.wrap_socket(request, server_side=True)           # the handshake runs here, in the request thread
        tls.settimeout(None)
        try:
            self.RequestHandlerClass(tls, client_address, self)
        finally:
            tls.close()

    def handle_error(self, request, client_address):
        e = sys.exc_info()[1]
        if isinstance(e, (ConnectionError, TimeoutError, ssl.SSLError)) or (isinstance(e, OSError) and e.errno in QUIET_ERRNO):
            if VERBOSE:
                print(f'{client_address[0]} - connection dropped: {type(e).__name__} {e}')
            return
        super().handle_error(request, client_address)


def bind(port, explicit):
    """the server on port, or on the next free port when the default one is taken"""
    for p in range(port, port + (1 if explicit else 20)):
        try:
            return Server(('0.0.0.0', p), Handler), p
        except OSError as e:
            if e.errno != errno.EADDRINUSE:
                raise
            print(f'port {p} is in use (another server still running? e.g. the Tsugite Joint FEA serve.py)'
                  + ('' if explicit else ', trying the next one'))
    sys.exit('no free port: stop the other server (Ctrl-C in its terminal) or pass --port')


def lan_addresses():
    addrs = set()
    try:
        for info in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET):
            ip = info[4][0]
            if not ip.startswith('127.'):
                addrs.add(ip)
    except socket.gaierror:
        pass
    try:                                                            # the interface that routes to the internet
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM); s.connect(('10.255.255.255', 1)); addrs.add(s.getsockname()[0]); s.close()
    except OSError:
        pass
    return sorted(addrs)


def lan_selftest(scheme, port, ips):
    """fetch the page through this Mac's own network address: if that fails, other devices cannot get in either"""
    time.sleep(0.5)
    ctx = ssl._create_unverified_context() if scheme == 'https' else None
    for ip in ips:
        url = f'{scheme}://{ip}:{port}/index.html'
        try:
            with urllib.request.urlopen(url, timeout=4, context=ctx) as r:
                ok = r.status == 200
        except Exception as e:                                          # noqa: BLE001 - any failure is the answer
            ok, why = False, f'{type(e).__name__}: {e}'
        if ok:
            print(f'  network check:  {ip} answers - if a phone still cannot connect, look at the Wi-Fi (guest network, '
                  f'client isolation) or the VPN; run with --verbose to see whether its requests arrive')
        else:
            print(f'  network check:  {ip} does NOT answer from this Mac ({why}).\n'
                  '    Likely: macOS blocks Python from the local network. Check, in this order:\n'
                  '    1. System Settings > Privacy & Security > Local Network: switch on the app this runs in (Terminal, iTerm, VS Code)\n'
                  f'    2. System Settings > Network > Firewall > Options: allow incoming connections for {sys.executable}\n'
                  '    3. a VPN (Cisco Secure Client) that blocks the local network: disconnect it at home\n'
                  '    then restart serve.py.')


def ensure_cert(path_pem):
    if os.path.exists(path_pem):
        return path_pem
    names = ','.join(['DNS:localhost', 'DNS:' + socket.gethostname()] + ['IP:' + ip for ip in lan_addresses()] + ['IP:127.0.0.1'])
    cmd = ['openssl', 'req', '-x509', '-newkey', 'rsa:2048', '-nodes', '-keyout', path_pem, '-out', path_pem, '-days', '825',
           '-subj', f'/CN={APP}', '-addext', 'subjectAltName=' + names]
    print('creating a self-signed certificate:', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    return path_pem


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--port', type=int, default=None)
    ap.add_argument('--https', action='store_true')
    ap.add_argument('--verbose', action='store_true', help='log every request and dropped connection')
    args = ap.parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    os.chdir(HERE)
    httpd, port = bind(args.port or (8443 if args.https else 8080), explicit=args.port is not None)
    scheme = 'http'
    if args.https:
        pem = ensure_cert(os.path.join(HERE, PEM))
        ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER); ctx.load_cert_chain(pem)
        httpd.tls = ctx; scheme = 'https'                            # TLS per connection, in its thread (Server.finish_request)
    print(f'{APP} — serving {HERE}')
    print(f'  this machine:   {scheme}://localhost:{port}/   <- open this one here (WebGPU needs localhost or https)')
    for ip in lan_addresses():
        print(f'  other devices:  {scheme}://{ip}:{port}/   (same Wi-Fi{"" if args.https else "; CPU solver only over plain http: use --https for WebGPU"})')
    print('Ctrl-C stops the server.')
    ips = lan_addresses()
    if ips:
        threading.Thread(target=lan_selftest, args=(scheme, port, ips), daemon=True).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print('\nstopped')


if __name__ == '__main__':
    sys.exit(main())
