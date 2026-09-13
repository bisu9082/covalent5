#!/usr/bin/env python3
"""verify_dois.py - resolve every DOI in references.bib against Crossref.

For each bibliography entry the stored DOI is resolved through the Crossref REST
API and the response is compared with the stored record on title (token-level
similarity) and publication year. Entries that fail to resolve are reported as
FAIL; entries that resolve to a record disagreeing with the stored metadata are
reported as WARN. Exit status is non-zero if any entry fails.

Usage
-----
    python verify_dois.py [references.bib]
"""
import json, re, subprocess, sys, time, difflib
UA = 'covalent5-verify/1.0 (mailto:bisu9082@gmail.com)'

def cr(doi):
    r = subprocess.run(['curl','-sS','--max-time','25','-H',f'User-Agent: {UA}',
                        f'https://api.crossref.org/works/{doi}'], capture_output=True, text=True)
    try:
        j = json.loads(r.stdout)
        return j['message'] if j.get('status') == 'ok' else None
    except Exception:
        return None

def norm(s):
    return re.sub(r'[^a-z0-9]', '', (s or '').lower())

def main(path='references.bib'):
    txt = open(path, encoding='utf-8').read()
    entries = [e for e in re.split(r'\n(?=@)', txt) if e.strip().startswith('@')]
    ok = warn = bad = 0
    print(f"{'key':<22}{'status':<9}DOI / note")
    print('-'*100)
    for e in entries:
        key = re.match(r'@\w+\{([^,]+),', e.strip()).group(1)
        dm = re.search(r'\bdoi\s*=\s*[{"]([^}"]+)', e, re.I)
        tm = re.search(r'\btitle\s*=\s*\{+([^}]+)', e, re.I)
        ym = re.search(r'\byear\s*=\s*[{"]?(\d{4})', e, re.I)
        if not dm:
            print(f"{key:<22}{'NO-DOI':<9}"); bad += 1; continue
        d = dm.group(1).strip()
        m = cr(d)
        time.sleep(0.4)
        if not m:
            print(f"{key:<22}{'FAIL':<9}{d}  <- not resolved by Crossref"); bad += 1; continue
        notes = []
        if tm:
            r = difflib.SequenceMatcher(None, norm(tm.group(1)), norm(m['title'][0])).ratio()
            if r < 0.80:
                notes.append(f"title mismatch ({r:.2f}): Crossref=\"{m['title'][0][:60]}\"")
        if ym:
            cy = m.get('published', {}).get('date-parts', [[None]])[0][0]
            if cy and abs(int(ym.group(1)) - cy) > 1:
                notes.append(f"year {ym.group(1)} vs Crossref {cy}")
        if notes:
            print(f"{key:<22}{'WARN':<9}{d}"); [print(' '*31 + n) for n in notes]; warn += 1
        else:
            print(f"{key:<22}{'ok':<9}{d}"); ok += 1
    print('-'*100)
    print(f"total {len(entries)}  |  ok {ok}  |  WARN {warn}  |  FAIL {bad}")
    return 0 if bad == 0 else 1

if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
