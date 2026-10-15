"""index.html (a complete document for local use) → artifact.html (body content only; the artifact host wraps it)."""
import re, os
here = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(here, 'index.html'), encoding='utf8').read()
s = re.sub(r'<!doctype html>\s*', '', s, flags=re.I)
s = re.sub(r'</?html[^>]*>\s*', '', s)
s = re.sub(r'<head>\s*', '', s); s = re.sub(r'</head>\s*', '', s)
s = re.sub(r'<body>\s*', '', s); s = re.sub(r'</body>\s*', '', s)
s = re.sub(r'<meta charset="utf-8">\s*', '', s)
s = re.sub(r'<meta name="viewport"[^>]*>\s*', '', s)
open(os.path.join(here, 'artifact.html'), 'w', encoding='utf8').write(s)
print('artifact.html', len(s), 'bytes; starts with:', s[:60].replace('\n', ' '))
